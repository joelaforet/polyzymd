"""Tests for the sasa and residue_sasa functions and ``polyzymd analyze sasa``.

The unit tests build MDAnalysis universes in memory with explicit elements.
The study tests write OpenMM run directories of a twelve-atom cluster in four
residues of three atoms, named ``C0`` to ``C11`` so the loader reads them as
carbon, and jittered at random from frame to frame.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import SASA_PROBE_RADIUS_NM, residue_sasa, sasa
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
md = pytest.importorskip("mdtraj")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

RESINDEX = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
EQUILIBRATION = "0ns"


def _cluster(seed: int = 0) -> np.ndarray:
    """Twelve atoms packed tightly enough that every residue occludes its neighbours."""
    return np.random.default_rng(seed).normal(scale=2.0, size=(12, 3))


def _universe(coordinates, elements, resindex=None) -> "mda.Universe":
    """An in-memory universe with the given elements and one frame per coordinate set."""
    coordinates = np.asarray(coordinates, dtype=np.float32)
    if coordinates.ndim == 2:
        coordinates = coordinates[np.newaxis]
    n_atoms = coordinates.shape[1]
    resindex = list(resindex if resindex is not None else range(n_atoms))
    universe = mda.Universe.empty(
        n_atoms, n_residues=max(resindex) + 1, atom_resindex=resindex, trajectory=True
    )
    universe.add_TopologyAttr("names", [f"A{i}" for i in range(n_atoms)])
    universe.add_TopologyAttr("resnames", ["ALA"] * (max(resindex) + 1))
    universe.add_TopologyAttr("resids", list(range(1, max(resindex) + 2)))
    universe.add_TopologyAttr("elements", list(elements))
    universe.load_new(coordinates, format="MEMORY")
    return universe


def _single_frame_atom_sasa(positions, elements) -> np.ndarray:
    """Per-atom SASA in Å² of one frame, from a fresh one-frame mdtraj call."""
    topology = md.Topology()
    residue = topology.add_residue("ALA", topology.add_chain())
    for index, symbol in enumerate(elements):
        topology.add_atom(f"A{index}", md.element.get_by_symbol(symbol), residue)
    xyz = np.asarray(positions, dtype=np.float32)[np.newaxis] / 10.0
    atom_nm2 = md.shrake_rupley(md.Trajectory(xyz=xyz, topology=topology), mode="atom")
    return np.asarray(atom_nm2, dtype=np.float64)[0] * 100.0


# ---------------------------------------------------------------------------
# The functions
# ---------------------------------------------------------------------------


def test_isolated_atom_has_the_full_sphere_area() -> None:
    """Every sphere point of an atom far from the others is accessible."""
    from mdtraj.geometry.sasa import _ATOMIC_RADII

    universe = _universe([[0, 0, 0], [100, 0, 0], [0, 100, 0]], ["C", "N", "O"])
    for atom, symbol in zip(universe.atoms, ("C", "N", "O"), strict=True):
        radius_a = 10.0 * (_ATOMIC_RADII[symbol] + SASA_PROBE_RADIUS_NM)
        expected = 4.0 * math.pi * radius_a**2
        assert sasa(universe.atoms[[atom.index]], universe.atoms) == pytest.approx(
            expected, rel=1e-6
        )
    assert sasa(universe.atoms, universe.atoms) == pytest.approx(
        sum(
            4.0 * math.pi * (10.0 * (_ATOMIC_RADII[s] + SASA_PROBE_RADIUS_NM)) ** 2
            for s in ("C", "N", "O")
        ),
        rel=1e-6,
    )


def test_overlapping_atoms_and_context_atoms_occlude_the_target() -> None:
    """Two atoms 2 Å apart expose less than twice one, and a context atom hides area."""
    universe = _universe([[0, 0, 0], [2, 0, 0], [100, 0, 0]], ["C", "C", "C"])
    pair, lone = universe.atoms[[0, 1]], universe.atoms[[2]]
    first = universe.atoms[[0]]

    assert sasa(pair, universe.atoms) < 2 * sasa(lone, universe.atoms)
    assert sasa(first, first) == pytest.approx(sasa(lone, universe.atoms))
    assert sasa(first, first) > sasa(first, pair)


def test_residue_sasa_computes_every_frame_alone() -> None:
    """Identical frames give exactly the single-frame per-residue SASA.

    MDTraj 1.11.1 gives a frame that follows another in the same
    ``shrake_rupley`` call (in the same OpenMP thread) about 0.1 percent too
    much area, so with frames batched this average would come out larger
    than the single-frame value for this cluster.
    """
    coordinates = _cluster()
    elements = ["C", "N", "O"] * 4
    universe = _universe(np.repeat(coordinates[np.newaxis], 48, axis=0), elements, RESINDEX)
    single = np.bincount(RESINDEX, weights=_single_frame_atom_sasa(coordinates, elements))

    result = residue_sasa(universe.atoms, universe.atoms, np.arange(48))

    assert result == pytest.approx(single, rel=1e-12, abs=0.0)
    assert residue_sasa(universe.atoms, universe.atoms, [0]) == pytest.approx(single, rel=1e-12)


def test_residue_sasa_follows_the_target_residues_in_a_larger_context() -> None:
    """Only target residues are reported, in order, and context atoms occlude them."""
    universe = _universe(_cluster(), ["C"] * 12, RESINDEX)
    target = universe.select_atoms("resid 2 3")

    alone = residue_sasa(target, target, [0])
    crowded = residue_sasa(target, universe.atoms, [0])

    assert alone.shape == crowded.shape == (2,)
    assert np.all(crowded <= alone) and crowded.sum() < alone.sum()
    assert crowded.sum() == pytest.approx(sasa(target, universe.atoms))


def test_same_context_with_another_target_is_measured_for_that_target() -> None:
    """The cached topology of a context is not reused with another target's atoms."""
    universe = _universe(_cluster(), ["C"] * 12, RESINDEX)
    first, second = universe.select_atoms("resid 1"), universe.select_atoms("resid 4")
    atom = _single_frame_atom_sasa(universe.atoms.positions, ["C"] * 12)

    assert sasa(first, universe.atoms) == pytest.approx(atom[:3].sum(), rel=1e-12)
    assert sasa(second, universe.atoms) == pytest.approx(atom[9:].sum(), rel=1e-12)


def test_target_outside_the_context_is_refused() -> None:
    universe = _universe(_cluster(), ["C"] * 12, RESINDEX)
    with pytest.raises(ProtocolError, match="must be a non-empty part of the context") as info:
        sasa(universe.select_atoms("resid 1 2"), universe.select_atoms("resid 2 3"))
    assert "contains every target atom" in info.value.hint
    with pytest.raises(ProtocolError, match="non-empty part"):
        sasa(universe.atoms[[]], universe.atoms)


def test_unknown_element_is_refused() -> None:
    universe = _universe([[0, 0, 0], [5, 0, 0]], ["C", "Qq"])
    with pytest.raises(ProtocolError, match=r"atom 1 \(A1\) has element 'Qq'"):
        sasa(universe.atoms, universe.atoms)


# ---------------------------------------------------------------------------
# The study API and polyzymd analyze sasa
# ---------------------------------------------------------------------------


def _frames(seed: int, jitter: float, n_frames: int = 24) -> np.ndarray:
    rng = np.random.default_rng(seed)
    base = _cluster()
    return np.array([base + rng.normal(scale=jitter, size=base.shape) for _ in range(n_frames)])


@pytest.fixture()
def coordinates() -> dict[tuple[str, int], np.ndarray]:
    return {
        (label, replicate): _frames(10 * replicate + len(label), jitter)
        for label, jitter in (("A", 0.1), ("B", 0.4))
        for replicate in (1, 2, 3)
    }


@pytest.fixture()
def configs(tmp_path: Path, coordinates) -> dict[str, Path]:
    """Two conditions of three replicates of the jittered cluster."""
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_frames(config, replicate, coordinates[(label, replicate)], RESINDEX)
        paths[label] = config
    return paths


def test_timeseries_values_equal_a_fresh_single_frame_mdtraj_call(
    configs, coordinates, tmp_path
) -> None:
    """Each stored frame is the SASA of that frame's coordinates computed on their own."""
    study = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
    series = study.timeseries(
        functions.sasa, pz.select("all"), pz.select("all"), unit="A^2", output_dir=tmp_path
    )

    for (label, replicate), frames in coordinates.items():
        stored = series.series[label][replicate - 1]
        expected = [
            _single_frame_atom_sasa(frames[k].astype(np.float32), ["C"] * 12).sum()
            for k in stored.frames
        ]
        assert stored.values == pytest.approx(expected, rel=1e-12)


def test_analyze_sasa_measures_the_target_alone_by_default(configs, tmp_path) -> None:
    report = analyze(
        "sasa",
        [configs["A"]],
        equilibration=EQUILIBRATION,
        settings={"target": "resid 1 2"},
        output_dir=tmp_path,
        plots=False,
    )

    assert (report.analysis, report.run, report.metric, report.unit) == (
        "sasa",
        "isolated",
        "mean_sasa",
        "A^2",
    )
    assert report.all_runs == ["isolated", "isolated_residues"]
    assert report.provenance.settings == {
        "target": "resid 1 2",
        "contexts": {"isolated": "resid 1 2"},
        "probe_radius_nm": 0.14,
        "n_sphere_points": 960,
    }
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    series = study.timeseries(
        functions.sasa,
        pz.select("resid 1 2"),
        pz.select("resid 1 2"),
        unit="A^2",
        output_dir=tmp_path / "direct",
    )
    means = [float(np.mean(item.values)) for item in series.series["A"]]
    (condition,) = report.conditions
    assert condition.replicate_values == pytest.approx(means, rel=1e-12)


def test_analyze_sasa_named_contexts_and_per_residue_results(configs, tmp_path) -> None:
    settings = {"target": "resid 1 2", "contexts": {"alone": "resid 1 2", "crowded": "all"}}
    options = {"equilibration": EQUILIBRATION, "settings": settings, "plots": False}
    both = [configs["A"], configs["B"]]

    alone = analyze("sasa", both, output_dir=tmp_path, **options)
    crowded = analyze("sasa", both, output_dir=tmp_path, run="crowded", **options)
    residues = analyze("sasa", both, output_dir=tmp_path, run="crowded_residues", **options)

    runs = ["alone", "alone_residues", "crowded", "crowded_residues"]
    assert alone.run == "alone" and alone.all_runs == crowded.all_runs == runs
    assert residues.provenance.settings["contexts"] == settings["contexts"]
    for label in ("A", "B"):
        (lone,) = [row for row in alone.conditions if row.label == label]
        (packed,) = [row for row in crowded.conditions if row.label == label]
        assert packed.mean < lone.mean
    assert residues.run == "crowded_residues"
    assert sorted((row.label, row.entry) for row in residues.conditions) == [
        (label, entry) for label in ("A", "B") for entry in ("1", "2")
    ]
    study = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
    direct = study.per_replicate(
        functions.residue_sasa,
        pz.select("resid 1 2"),
        pz.select("all"),
        unit="A^2",
        labels=lambda u: u.select_atoms("resid 1 2").residues.resids,
        output_dir=tmp_path / "direct",
    )
    for row in residues.conditions:
        position = int(row.entry) - 1
        expected = [float(values[position]) for values in direct.values[row.label]]
        assert row.replicate_values == pytest.approx(expected, rel=1e-12)


def test_analyze_sasa_refuses_an_unknown_run_and_a_bad_context(configs) -> None:
    settings = {"target": "resid 1 2", "contexts": {"alone": "resid 1 2"}}
    with pytest.raises(ProtocolError, match="no result named 'crowded'") as info:
        analyze(
            "sasa", [configs["A"]], equilibration=EQUILIBRATION, settings=settings, run="crowded"
        )
    assert "['alone', 'alone_residues']" in info.value.hint
    with pytest.raises(ProtocolError, match="contexts must map names to selections"):
        analyze("sasa", [configs["A"]], settings={"contexts": ["all"]}, plots=False)
    with pytest.raises(ProtocolError, match="non-empty part of the context"):
        analyze(
            "sasa",
            [configs["A"]],
            equilibration=EQUILIBRATION,
            settings={"target": "all", "contexts": {"part": "resid 1"}},
            plots=False,
        )


def test_cli_sasa_draws_the_documented_figures(configs, tmp_path) -> None:
    arguments = ["sasa", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--set", "target=resid 1 2", "--set", "contexts={crowded: all}"]
    arguments += ["--output-dir", str(tmp_path)]

    total = CliRunner().invoke(analyze_command, arguments)
    per_residue = CliRunner().invoke(analyze_command, [*arguments, "--run", "crowded_residues"])

    assert total.exit_code == 0, total.output
    assert total.stdout.startswith("# polyzymd analyze sasa  metric mean_sasa  unit A^2")
    assert per_residue.exit_code == 0, per_residue.output
    rows = [line for line in per_residue.stdout.split("\n") if line.startswith("A vs B  ")]
    assert rows[0].startswith("A vs B  labels 2  tested 2  family 2")
    assert {path.name for path in (tmp_path / "figures" / "sasa").iterdir()} == {
        "sasa_timeseries_crowded.png",
        "sasa_comparison_crowded.png",
        "sasa_distribution_crowded.png",
        "sasa_profile_crowded.png",
        "sasa_difference_crowded.png",
    }
