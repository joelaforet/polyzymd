"""Tests for the contract SASA plugin."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.sasa import SASA, SASAAnalysis, SASARun, SASASettings


def make_two_residue_universe(separation: float = 40.0, n_frames: int = 3) -> Any:
    """Two carbon atoms in separate alanine residues, far enough apart to be isolated.

    Parameters
    ----------
    separation : float, optional
        Distance between the two atoms in angstrom, by default 40.0.
    n_frames : int, optional
        Number of identical frames, by default 3.

    Returns
    -------
    MDAnalysis.Universe
        Universe backed by ``MemoryReader``, writable as a PDB file.
    """
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(2, n_residues=2, atom_resindex=[0, 1], trajectory=True)
    universe.add_TopologyAttr("names", ["CA", "CA"])
    universe.add_TopologyAttr("types", ["C", "C"])
    universe.add_TopologyAttr("elements", ["C", "C"])
    universe.add_TopologyAttr("resnames", ["ALA", "ALA"])
    universe.add_TopologyAttr("resids", [1, 2])
    universe.add_TopologyAttr("segids", ["A"])
    positions = np.asarray([[0.0, 0.0, 0.0], [separation, 0.0, 0.0]], dtype=np.float32)
    universe.load_new(np.stack([positions] * n_frames), format=MemoryReader)
    return universe


#: Area of one isolated MDTraj carbon sphere (radius 0.17 nm) with a 0.14 nm probe,
#: in square angstrom. Every sphere point is accessible, so Shrake-Rupley returns
#: 4 pi r^2 exactly.
ISOLATED_CARBON_A2 = 4.0 * np.pi * (0.17 + 0.14) ** 2 * 100.0

SETTINGS = SASASettings(
    runs=[SASARun(label="all", target_selection="all", context_selection="all")]
)


def make_cluster_universe(n_atoms: int = 30, n_frames: int = 40) -> Any:
    """One overlapping cluster of carbons, repeated as identical frames.

    Parameters
    ----------
    n_atoms : int, optional
        Atoms in the cluster, by default 30.
    n_frames : int, optional
        Number of identical copies of the one frame, by default 40.

    Returns
    -------
    MDAnalysis.Universe
        Universe whose every frame holds the same coordinates.
    """
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        n_atoms, n_residues=n_atoms, atom_resindex=list(range(n_atoms)), trajectory=True
    )
    universe.add_TopologyAttr("names", ["CA"] * n_atoms)
    universe.add_TopologyAttr("types", ["C"] * n_atoms)
    universe.add_TopologyAttr("elements", ["C"] * n_atoms)
    universe.add_TopologyAttr("resnames", ["ALA"] * n_atoms)
    universe.add_TopologyAttr("resids", list(range(1, n_atoms + 1)))
    universe.add_TopologyAttr("segids", ["A"])
    frame = np.random.default_rng(0).normal(scale=3.5, size=(n_atoms, 3)).astype(np.float32)
    universe.load_new(np.stack([frame] * n_frames), format=MemoryReader)
    return universe


@pytest.fixture
def two_residues() -> Any:
    """Two isolated carbon atoms in separate alanine residues."""
    return make_two_residue_universe()


class TestObservables:
    """The plugin reports one total and one per-residue profile per context."""

    def test_isolated_spheres_match_the_analytic_area(self, two_residues: Any) -> None:
        """Two atoms that never touch have the analytic area of two spheres."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        observables = SASA().compute(
            two_residues, FrameSelection(start=0, stop=3, step=1), SETTINGS
        )
        total = next(obs for obs in observables if obs.name == "sasa_all")

        assert total.unit == "A^2"
        assert total.kind == "mean_of_timeseries"
        assert total.metadata["chunk_size"] == 100
        assert total.metadata["probe_radius_nm"] == 0.14
        assert total.metadata["n_sphere_points"] == 960
        np.testing.assert_allclose(total.values, [2.0 * ISOLATED_CARBON_A2] * 3, rtol=1e-6)

    def test_relative_profile_divides_by_the_tien_maximum(self, two_residues: Any) -> None:
        """Relative area is the residue area over the Tien et al. 2013 maximum."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        observables = SASA().compute(
            two_residues, FrameSelection(start=0, stop=3, step=1), SETTINGS
        )
        profile = next(obs for obs in observables if obs.name == "relative_sasa_all")

        assert profile.kind == "profile"
        assert profile.unit == "fraction"
        assert profile.index == [0.0, 1.0]
        assert profile.index_label == "residue index"
        assert profile.metadata["residue_labels"] == ["A:1:ALA", "A:2:ALA"]
        np.testing.assert_allclose(profile.values, [ISOLATED_CARBON_A2 / 121.0] * 2, rtol=1e-6)

    def test_protonation_variant_uses_its_parent_residue(self) -> None:
        """A force field variant name such as HIE takes the HIS maximum area."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        universe = make_two_residue_universe()
        universe.residues.resnames = ["HIE", "HIS"]

        observables = SASA().compute(universe, FrameSelection(start=0, stop=3, step=1), SETTINGS)
        profile = next(obs for obs in observables if obs.name == "relative_sasa_all")

        np.testing.assert_allclose(profile.values, [ISOLATED_CARBON_A2 / 216.0] * 2, rtol=1e-6)


class TestLifecycle:
    """The framework aggregates the contract observables without plugin help."""

    def test_aggregate_reports_uncertainty_over_replicates(
        self, two_residues: Any, run_contract_analysis: Any
    ) -> None:
        """Three identical replicates give a zero SEM on the total area."""
        artifact = run_contract_analysis(SASAAnalysis, SETTINGS, two_residues)
        aggregates = {
            payload["name"]: ObservableAggregate.model_validate(payload)
            for payload in artifact.payload["observables"]
        }

        total = aggregates["sasa_all"]
        assert total.n_replicates == 3
        assert total.mean == pytest.approx(2.0 * ISOLATED_CARBON_A2, rel=1e-6)
        assert total.sem == pytest.approx(0.0)
        profile = aggregates["relative_sasa_all"]
        assert profile.index == [0.0, 1.0]
        assert profile.index_label == "residue index"
        assert profile.metadata["chunk_size"] == 100

    def test_analysis_keeps_the_expensive_resource_hints(self) -> None:
        """The generated class carries the hints the orchestrator submits with."""
        assert SASAAnalysis.execution_cost_hint == "high"
        assert SASAAnalysis.slurm_resource_hint.mem == "8G"


class TestChunkSize:
    """chunk_size changes the numbers, so it is pinned and recorded."""

    def test_chunking_changes_the_totals_by_about_a_tenth_of_a_percent(self) -> None:
        """MDTraj gives one frame different areas depending on the array it is in.

        ``mdtraj.shrake_rupley`` takes a different code path once the frame
        array is long enough, so identical frames come back with slightly
        different areas. The effect is real and small, and it is pinned here so
        a future MDTraj release that removes it is noticed rather than quietly
        changing every stored SASA number. If this test fails because the two
        series are now equal, delete the warnings about chunk_size in the
        module docstring and in the reference page.
        """
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        universe = make_cluster_universe()
        frames = FrameSelection(start=0, stop=40, step=1)

        def totals(chunk_size: int) -> np.ndarray:
            settings = SASASettings(
                runs=[SASARun(label="all", target_selection="all")], chunk_size=chunk_size
            )
            observables = SASA().compute(universe, frames, settings)
            return np.asarray(next(o for o in observables if o.name == "sasa_all").values)

        one_at_a_time, all_at_once = totals(1), totals(40)

        assert not np.array_equal(one_at_a_time, all_at_once)
        spread = np.max(np.abs(all_at_once - one_at_a_time)) / np.mean(one_at_a_time)
        assert spread < 2e-3

    def test_one_frame_per_chunk_is_self_consistent(self) -> None:
        """With one frame per call, identical frames give identical areas."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        settings = SASASettings(runs=[SASARun(label="all", target_selection="all")], chunk_size=1)
        observables = SASA().compute(
            make_cluster_universe(), FrameSelection(start=0, stop=40, step=1), settings
        )
        values = np.asarray(next(o for o in observables if o.name == "sasa_all").values)

        assert len(np.unique(values)) == 1


class TestInvalidInput:
    """An input the method does not apply to raises rather than returning zero."""

    def test_empty_selection_raises(self, two_residues: Any) -> None:
        """A selection matching no atoms is an error, not an area of zero."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        settings = SASASettings(runs=[SASARun(label="x", target_selection="resname NOPE")])

        with pytest.raises(SelectionError, match="matched no atoms"):
            SASA().compute(two_residues, FrameSelection(start=0, stop=3, step=1), settings)

    def test_target_outside_context_raises(self, two_residues: Any) -> None:
        """A target the context does not contain would lose its own neighbours."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        settings = SASASettings(
            runs=[SASARun(label="x", target_selection="all", context_selection="resid 1")]
        )

        with pytest.raises(ReplicateError, match="outside the context selection"):
            SASA().compute(two_residues, FrameSelection(start=0, stop=3, step=1), settings)

    def test_non_standard_residue_has_no_relative_area(self) -> None:
        """A residue outside the Tien table has no defined relative area."""
        from polyzymd.analyses.mda.frame_selection import FrameSelection

        universe = make_two_residue_universe()
        universe.residues.resnames = ["SBM", "SBM"]

        with pytest.raises(ReplicateError, match="Tien"):
            SASA().compute(universe, FrameSelection(start=0, stop=3, step=1), SETTINGS)

    def test_duplicate_labels_rejected(self) -> None:
        """Two contexts with one label would collide on the observable name."""
        with pytest.raises(ValueError, match="labels must be unique"):
            SASASettings(
                runs=[
                    SASARun(label="a", target_selection="protein"),
                    SASARun(label="a", target_selection="protein"),
                ]
            )

    def test_per_run_stride_warns_where_the_user_will_see_it(self) -> None:
        """A per-context stride changes how much data is analysed, so it warns loudly."""
        with pytest.warns(UserWarning, match="ignored since"):
            SASARun(label="a", target_selection="protein", stride=5)
