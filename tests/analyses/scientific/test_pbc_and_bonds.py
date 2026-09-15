"""Known-answer tests for periodic boundary handling and bond requirements.

These tests pin three behaviours. Fragment-based observables (Rg in fragment
mode, contacts chain identity) must fail loudly when the topology carries no
bonds instead of collapsing the selection into a single fragment. Pair
distances must use the minimum image convention against the box of the frame
they were measured in, which means no rigid-body alignment may run first.
Universe provenance must state which periodic boundary policy was applied.
"""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")

BOX_LENGTH = 50.0


def _two_chain_universe() -> Any:
    """Build a two-chain polymer universe with explicit bonds.

    Returns
    -------
    Any
        MDAnalysis universe with six atoms in two bonded chains.
    """

    universe = mda.Universe.empty(
        6,
        n_residues=2,
        atom_resindex=[0, 0, 0, 1, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["C1", "C2", "C3", "C1", "C2", "C3"])
    universe.add_TopologyAttr("type", ["C"] * 6)
    universe.add_TopologyAttr("resname", ["SBM", "SBM"])
    universe.add_TopologyAttr("resid", [1, 2])
    universe.add_TopologyAttr("segid", ["C"])
    universe.add_TopologyAttr("mass", [12.0] * 6)
    universe.add_bonds([(0, 1), (1, 2), (3, 4), (4, 5)])
    universe.atoms.positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [14.0, 0.0, 0.0],
            [18.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    )
    universe.dimensions = [BOX_LENGTH, BOX_LENGTH, BOX_LENGTH, 90.0, 90.0, 90.0]
    return universe


def _fragment_run(**overrides: Any) -> Any:
    """Build an Rg run in fragment mode.

    Parameters
    ----------
    **overrides : Any
        Field overrides for ``RgRunSettings``.

    Returns
    -------
    Any
        Configured ``RgRunSettings`` instance.
    """

    from polyzymd.analyses.rg import RgRunSettings

    fields: dict[str, Any] = {
        "label": "polymer",
        "selection": "all",
        "calculation_mode": "fragments",
        "save_fragment_distribution": False,
    }
    fields.update(overrides)
    return RgRunSettings(**fields)


def _run_rg(universe: Any, run: Any) -> dict[str, Any]:
    """Measure one Rg run over a whole universe, keyed by observable name.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe to measure.
    run : Any
        Configured Rg run settings.

    Returns
    -------
    dict[str, Any]
        Observables the plugin returned.
    """

    from polyzymd.analyses.mda.frame_selection import FrameSelection
    from polyzymd.analyses.rg import Rg, RgSettings

    frames = FrameSelection(start=0, stop=len(universe.trajectory), step=1)
    return {
        observable.name: observable
        for observable in Rg().compute(universe, frames, RgSettings(runs=[run]))
    }


def test_rg_fragment_mode_measures_each_bonded_chain() -> None:
    """Fragment mode should report one Rg per bonded chain."""

    universe = _two_chain_universe()

    observables = _run_rg(universe, _fragment_run())

    fragment_rg = np.asarray(observables["rg_polymer_fragments"].values)
    assert fragment_rg.shape == (2,)
    np.testing.assert_allclose(
        fragment_rg,
        [math.sqrt(2.0 / 3.0), math.sqrt(32.0 / 3.0)],
        rtol=1e-6,
    )
    whole_selection_rg = float(universe.atoms.radius_of_gyration())
    reduced = float(observables["rg_polymer"].values[0])
    assert not math.isclose(reduced, whole_selection_rg, rel_tol=1e-3)


def test_rg_fragment_mode_without_bonds_raises_typed_error() -> None:
    """Fragment mode must refuse a topology with no bonds."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    universe = _two_chain_universe()
    universe.del_TopologyAttr("bonds")

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        _run_rg(universe, _fragment_run())

    message = str(excinfo.value)
    assert "6 atoms" in message
    assert "guess bonds" in message
    assert "bond" in message.lower()


def test_rg_fragment_mode_fallback_is_opt_in() -> None:
    """The old single-fragment fallback stays reachable behind a flag."""

    universe = _two_chain_universe()
    universe.del_TopologyAttr("bonds")

    observables = _run_rg(universe, _fragment_run(allow_single_fragment_fallback=True))

    assert len(observables["rg_polymer_fragments"].values) == 1
    assert "fragment_fallback" in observables["rg_polymer"].metadata


def test_identify_polymer_chains_without_bonds_raises_typed_error() -> None:
    """Contacts chain identity must refuse a topology with no bonds."""

    from polyzymd.analyses.contacts import ContactsSettings, _polymer_chains
    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    universe = _two_chain_universe()
    universe.del_TopologyAttr("bonds")

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        _polymer_chains(universe.atoms, ContactsSettings())

    assert "6 atoms" in str(excinfo.value)


def test_identify_polymer_chains_uses_fragments_when_bonds_exist() -> None:
    """Chain identity should assign one chain index per bonded fragment."""

    from polyzymd.analyses.contacts import ContactsSettings, _polymer_chains

    universe = _two_chain_universe()

    assert _polymer_chains(universe.atoms, ContactsSettings()).tolist() == [0, 1]


def _partially_bonded_universe() -> Any:
    """Build a universe where only the protein carries bonds.

    This is what a PDB with usable CONECT records for standard residues and
    none for the polymer produces. MDAnalysis raises nothing here, so an
    unguarded fragment call returns one singleton fragment per polymer atom.

    Returns
    -------
    Any
        MDAnalysis universe with a bonded three-atom protein residue and an
        unbonded three-atom polymer residue.
    """

    universe = mda.Universe.empty(
        6,
        n_residues=2,
        atom_resindex=[0, 0, 0, 1, 1, 1],
        residue_segindex=[0, 1],
        n_segments=2,
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["N", "CA", "C", "C1", "C2", "C3"])
    universe.add_TopologyAttr("type", ["N", "C", "C", "C", "C", "C"])
    universe.add_TopologyAttr("resname", ["ALA", "SBM"])
    universe.add_TopologyAttr("resid", [1, 2])
    universe.add_TopologyAttr("segid", ["A", "C"])
    universe.add_TopologyAttr("chainID", ["A", "A", "A", "C", "C", "C"])
    universe.add_TopologyAttr("mass", [14.0, 12.0, 12.0, 12.0, 12.0, 12.0])
    universe.add_bonds([(0, 1), (1, 2)])
    universe.atoms.positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [14.0, 0.0, 0.0],
            [18.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    )
    universe.dimensions = [BOX_LENGTH, BOX_LENGTH, BOX_LENGTH, 90.0, 90.0, 90.0]
    return universe


def test_partially_bonded_topology_is_detected_as_unbonded_selection() -> None:
    """A selection with no bonds of its own fails even when the protein has bonds."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError
    from polyzymd.analyses.shared.topology import topology_bond_source

    universe = _partially_bonded_universe()

    # The universe as a whole does carry bonds, so the topology-level report
    # cannot catch this on its own.
    assert topology_bond_source(universe) == (True, "conect")
    assert [len(fragment) for fragment in universe.select_atoms("resname SBM").fragments] == [
        1,
        1,
        1,
    ]

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        _run_rg(universe, _fragment_run(selection="resname SBM"))

    assert "3 atoms" in str(excinfo.value)


def test_partially_bonded_topology_rejects_polymer_chain_identity() -> None:
    """Contacts chain identity fails on a polymer selection with no bonds."""

    from polyzymd.analyses.contacts import ContactsSettings, _polymer_chains
    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    universe = _partially_bonded_universe()

    with pytest.raises(TopologyBondsMissingError):
        _polymer_chains(universe.select_atoms("resname SBM"), ContactsSettings())


def test_partially_bonded_topology_allows_the_bonded_selection() -> None:
    """The bonded protein selection still resolves its fragment."""

    universe = _partially_bonded_universe()

    observables = _run_rg(universe, _fragment_run(selection="resname ALA"))

    assert len(observables["rg_polymer_fragments"].values) == 1


def _mostly_bonded_universe(n_unbonded: int) -> Any:
    """Build a 40-atom universe with a chosen number of unbonded atoms.

    Parameters
    ----------
    n_unbonded : int
        How many of the 40 atoms are left without bonds.

    Returns
    -------
    Any
        MDAnalysis universe whose first ``40 - n_unbonded`` atoms form one
        bonded chain.
    """

    n_atoms = 40
    n_bonded = n_atoms - n_unbonded
    universe = mda.Universe.empty(
        n_atoms,
        n_residues=1,
        atom_resindex=[0] * n_atoms,
        residue_segindex=[0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", [f"C{index}" for index in range(n_atoms)])
    universe.add_TopologyAttr("type", ["C"] * n_atoms)
    universe.add_TopologyAttr("resname", ["SBM"])
    universe.add_TopologyAttr("resid", [1])
    universe.add_TopologyAttr("segid", ["C"])
    universe.add_TopologyAttr("mass", [12.0] * n_atoms)
    universe.add_bonds([(index, index + 1) for index in range(n_bonded - 1)])
    positions = np.zeros((n_atoms, 3), dtype=np.float32)
    positions[:, 0] = np.arange(n_atoms, dtype=np.float32)
    universe.atoms.positions = positions
    universe.dimensions = [BOX_LENGTH, BOX_LENGTH, BOX_LENGTH, 90.0, 90.0, 90.0]
    return universe


def test_a_few_unbonded_atoms_are_tolerated() -> None:
    """One stray atom in forty is below the singleton limit."""

    from polyzymd.analyses.shared.topology import require_topology_bonds

    fragments, fallback = require_topology_bonds(
        _mostly_bonded_universe(1).atoms, context="Rg run 'polymer' in fragment mode"
    )

    assert fallback is None
    assert len(fragments) == 2


def test_mostly_unbonded_selection_is_rejected_with_counts() -> None:
    """A selection where singletons pass the limit fails and names the counts."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError
    from polyzymd.analyses.shared.topology import require_topology_bonds

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        require_topology_bonds(
            _mostly_bonded_universe(4).atoms,
            context="Rg run 'polymer' in fragment mode",
        )

    message = str(excinfo.value)
    assert "4 of the 40 selected atoms are single-atom fragments" in message
    assert "10 percent" in message


def test_partially_bonded_selection_is_rejected_by_rg_fragment_mode() -> None:
    """A protein-plus-polymer selection with unbonded polymer fails, not averages to zero."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    universe = _partially_bonded_universe()

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        _run_rg(universe, _fragment_run(selection="all"))

    assert "3 of the 6 selected atoms are single-atom fragments" in str(excinfo.value)


def _pair_universe() -> Any:
    """Build a two-atom universe whose pair spans the periodic boundary.

    Returns
    -------
    Any
        MDAnalysis universe with a 50 Angstrom cubic box.
    """

    universe = mda.Universe.empty(
        2,
        n_residues=2,
        n_segments=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 1],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["CA", "CA"])
    universe.add_TopologyAttr("type", ["C", "C"])
    universe.add_TopologyAttr("resname", ["ALA", "ALA"])
    universe.add_TopologyAttr("resid", [1, 2])
    universe.add_TopologyAttr("segid", ["A", "B"])
    universe.add_TopologyAttr("mass", [12.0, 12.0])
    universe.atoms.positions = np.array(
        [[1.0, 0.0, 0.0], [49.0, 0.0, 0.0]],
        dtype=np.float32,
    )
    universe.dimensions = [BOX_LENGTH, BOX_LENGTH, BOX_LENGTH, 90.0, 90.0, 90.0]
    return universe


def _pair(label: str = "pair") -> list[Any]:
    """Build one pair selection for the two-atom universe."""

    from polyzymd.analyses.mda.pair_distance import PairSelection

    return [PairSelection(label=label, selection_a="index 0", selection_b="index 1")]


def _all_frames() -> Any:
    """Frame selection covering the whole trajectory."""

    from polyzymd.analyses.mda import FrameSelection

    return FrameSelection(start=0, stop=None, step=1, timestep_ps=1.0)


@pytest.mark.parametrize(("use_pbc", "expected"), [(True, 2.0), (False, 48.0)])
def test_pair_distance_minimum_image_matches_known_answer(use_pbc: bool, expected: float) -> None:
    """Minimum image folds the 48 Angstrom separation to 2 Angstrom."""

    from polyzymd.analyses.mda.pair_distance import pair_distance_matrix

    universe = _pair_universe()
    matrix = pair_distance_matrix(universe, _all_frames(), _pair(), use_pbc=use_pbc)

    np.testing.assert_allclose(matrix[0][0], expected, atol=1e-5)


def _rotate_about_z(universe: Any, degrees: float) -> None:
    """Rotate every coordinate about z while leaving the box unchanged.

    This reproduces what ``AlignTraj(..., in_memory=True)`` does to a universe.

    Parameters
    ----------
    universe : Any
        Universe whose coordinates are rotated in place.
    degrees : float
        Rotation angle in degrees.
    """

    angle = math.radians(degrees)
    rotation = np.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    universe.atoms.positions = universe.atoms.positions @ rotation.T


def test_rotating_coordinates_breaks_minimum_image_distances() -> None:
    """A rigid rotation with an unrotated box changes the folded distance."""

    from polyzymd.analyses.mda.pair_distance import pair_distance_matrix

    universe = _pair_universe()
    _rotate_about_z(universe, 45.0)
    matrix = pair_distance_matrix(universe, _all_frames(), _pair(), use_pbc=True)

    assert not math.isclose(float(matrix[0][0]), 2.0, abs_tol=1e-3)


def test_distances_plugin_keeps_minimum_image_and_refuses_alignment() -> None:
    """The distances plugin measures raw coordinates and deprecates alignment."""

    from polyzymd.analyses.distances import Distances, DistancesSettings

    with pytest.warns(UserWarning, match="align_trajectory"):
        settings = DistancesSettings(
            pairs=[{"label": "pair", "selection_a": "index 0", "selection_b": "index 1"}],
            use_pbc=True,
            align_trajectory=True,
        )
    observables = Distances().compute(_pair_universe(), _all_frames(), settings)

    np.testing.assert_allclose(observables[0].values[0], 2.0, atol=1e-5)


def test_triad_plugin_keeps_minimum_image() -> None:
    """The catalytic triad plugin measures raw coordinates too."""

    from polyzymd.analyses.catalytic_triad import CatalyticTriad, CatalyticTriadSettings

    settings = CatalyticTriadSettings(
        pairs=[{"label": "pair", "selection_a": "index 0", "selection_b": "index 1"}]
    )
    observables = CatalyticTriad().compute(_pair_universe(), _all_frames(), settings)

    np.testing.assert_allclose(observables[0].values[0], 2.0, atol=1e-5)


class _ProvenanceLoader:
    """Trajectory loader stub that records the requested PBC policy."""

    def __init__(self, info: Any, universe: Any) -> None:
        """Store the metadata and universe this stub returns.

        Parameters
        ----------
        info : Any
            Trajectory metadata returned by ``get_trajectory_info``.
        universe : Any
            Universe returned by ``load_universe``.
        """

        self.info = info
        self.universe = universe
        self.policies: list[str] = []

    def get_trajectory_info(self, replicate: int) -> Any:
        """Return the stored trajectory metadata."""

        return self.info

    def load_universe(
        self,
        replicate: int,
        cache: bool = True,
        *,
        pbc_policy: str = "as_is",
    ) -> Any:
        """Return the stored universe and record the requested policy."""

        self.policies.append(pbc_policy)
        return self.universe


def _trajectory_info(tmp_path: Path, trajectory_name: str = "prod_centered.xtc") -> Any:
    """Build trajectory metadata backed by real files.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory holding the fake inputs.
    trajectory_name : str, optional
        Trajectory filename, by default "prod_centered.xtc".

    Returns
    -------
    Any
        Populated ``TrajectoryInfo``.
    """

    from polyzymd.analyses.shared.loader import TrajectoryInfo

    working_dir = tmp_path / "run_1"
    working_dir.mkdir(parents=True, exist_ok=True)
    topology = working_dir / "solvated_system.pdb"
    topology.write_bytes(b"ATOM")
    trajectory = working_dir / trajectory_name
    trajectory.write_bytes(b"XTC")
    return TrajectoryInfo(
        topology_file=topology,
        trajectory_files=[trajectory],
        n_segments=1,
        working_directory=working_dir,
        replicate=1,
        topology_format="pdb",
        trajectory_format="xtc",
    )


def test_universe_provenance_records_pbc_policy_and_bond_source(tmp_path: Path) -> None:
    """Provenance must state the PBC policy, bond source, and trajectory variant."""

    from polyzymd.analyses.mda import UniverseProvider

    universe = _two_chain_universe()
    loader = _ProvenanceLoader(_trajectory_info(tmp_path), universe)
    provider = UniverseProvider(config=SimpleNamespace(engine="gromacs"), loader=loader)

    provider.load_universe(1)
    provenance = provider.provenance_for(1)

    assert provenance.pbc_policy == "as_is"
    assert provenance.topology_has_bonds is True
    assert provenance.bond_source == "conect"
    assert provenance.trajectory_variant == "centered"
    payload = provenance.as_dict()
    assert payload["pbc_policy"] == "as_is"
    assert payload["topology_has_bonds"] is True
    assert payload["bond_source"] == "conect"
    assert payload["trajectory_variant"] == "centered"


def test_universe_provider_forwards_make_whole_policy(tmp_path: Path) -> None:
    """A make_whole request reaches the loader and is recorded in provenance."""

    from polyzymd.analyses.mda import UniverseProvider

    universe = _two_chain_universe()
    loader = _ProvenanceLoader(_trajectory_info(tmp_path, "prod.xtc"), universe)
    provider = UniverseProvider(
        config=SimpleNamespace(engine="gromacs"),
        loader=loader,
        pbc_policy="make_whole",
    )

    provider.load_universe(1)
    provenance = provider.provenance_for(1)

    assert loader.policies == ["make_whole"]
    assert provenance.pbc_policy == "make_whole"
    assert provenance.trajectory_variant == "raw"


def test_provenance_refresh_keeps_bond_facts(tmp_path: Path) -> None:
    """Rediscovering input files must not forget what loading established."""

    from polyzymd.analyses.mda import UniverseProvider

    universe = _two_chain_universe()
    loader = _ProvenanceLoader(_trajectory_info(tmp_path), universe)
    provider = UniverseProvider(config=SimpleNamespace(engine="gromacs"), loader=loader)

    provider.load_universe(1)
    refreshed = provider.provenance_for(1, refresh=True)

    assert refreshed.topology_has_bonds is True
    assert refreshed.bond_source == "conect"


def test_every_condition_failing_the_same_way_is_reported(tmp_path: Path) -> None:
    """A typed error shared by every condition reaches the raised message."""

    from polyzymd.analyses._framework.lifecycle import _no_conditions_message
    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    error = TopologyBondsMissingError(
        context="contacts polymer chain detection",
        n_atoms=512000,
        topology=tmp_path / "solvated_system.pdb",
    )
    message = _no_conditions_message("contacts", [("no_polymer", error), ("sbma", error)])

    assert "TopologyBondsMissingError" in message
    assert "guess bonds" in message
    assert "512000 atoms" in message

    mixed = _no_conditions_message("contacts", [("a", error), ("b", ValueError("other"))])
    assert mixed == "contacts: no conditions succeeded analysis."


def test_loader_make_whole_requires_bonds(tmp_path: Path) -> None:
    """``make_whole`` on a bond-free topology raises the typed error."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError
    from polyzymd.analyses.shared.loader import apply_pbc_policy

    universe = _two_chain_universe()
    universe.del_TopologyAttr("bonds")

    with pytest.raises(TopologyBondsMissingError):
        apply_pbc_policy(universe, "make_whole", topology=tmp_path / "solvated_system.pdb")


def test_loader_make_whole_checks_the_selection_it_unwraps(tmp_path: Path) -> None:
    """``make_whole`` fails when most of its own selection has no bonds."""

    from polyzymd.analyses.exceptions import TopologyBondsMissingError
    from polyzymd.analyses.shared.loader import apply_pbc_policy
    from polyzymd.analyses.shared.topology import topology_bond_source

    universe = _partially_bonded_universe()

    # The topology-level report says bonds exist, which is why the check has to
    # be made against the atoms that are about to be unwrapped.
    assert topology_bond_source(universe)[0] is True

    with pytest.raises(TopologyBondsMissingError) as excinfo:
        apply_pbc_policy(universe, "make_whole", topology=tmp_path / "solvated_system.pdb")

    assert "single-atom fragments" in str(excinfo.value)


def test_loader_make_whole_unwraps_split_molecule(tmp_path: Path) -> None:
    """``make_whole`` rejoins a chain that straddles the periodic boundary."""

    from polyzymd.analyses.shared.loader import apply_pbc_policy

    universe = _two_chain_universe()
    positions = universe.atoms.positions.copy()
    positions[2] = [BOX_LENGTH - 1.0, 0.0, 0.0]
    universe.atoms.positions = positions
    split_rg = float(universe.select_atoms("index 0 1 2").radius_of_gyration())

    apply_pbc_policy(universe, "make_whole", topology=tmp_path / "solvated_system.pdb")
    universe.trajectory[0]

    whole_rg = float(universe.select_atoms("index 0 1 2").radius_of_gyration())
    assert whole_rg < split_rg
