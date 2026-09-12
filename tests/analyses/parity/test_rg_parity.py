"""Parity of the ported rg plugin against the pre-port implementation.

The reference in ``rg_reference.json`` was frozen from ``rg/_mda.py`` as it
stood before the port, on frames 500 to 600 of the 50:50 SBMA-EGMA run. Both
implementations call the same ``AtomGroup.radius_of_gyration`` on the same
frames, so the values must agree to within floating point noise.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.rg import Rg, RgRunSettings, RgSettings

REFERENCE_PATH = Path(__file__).with_name("rg_reference.json")
CONTROL_TOPOLOGY = Path(
    "/home/joelaforet/Shirts-Lab-Linux/polyzymd-realdata"
    "/LipA_ResorufinButyrate_noPoly_1000ns_363K_run1/solvated_system.pdb"
)
CONTROL_TRAJECTORY = CONTROL_TOPOLOGY.parent / "production_0" / "production_0_trajectory.dcd"
ATOL = 1e-8


@pytest.fixture(scope="module")
def reference() -> dict[str, Any]:
    """Frozen values from the pre-port implementation."""
    return json.loads(REFERENCE_PATH.read_text())


def _universe(topology: str | Path, trajectory: str | Path) -> Any:
    """Load a real trajectory, skipping the test when it is not on this machine."""
    if not Path(topology).exists() or not Path(trajectory).exists():
        pytest.skip(f"real trajectory not present at {topology}")
    mda = pytest.importorskip("MDAnalysis")
    return mda.Universe(str(topology), str(trajectory))


@pytest.fixture(scope="module")
def observables(reference: dict[str, Any]) -> dict[str, Any]:
    """Observables the ported plugin computes on the frozen window."""
    universe = _universe(reference["topology"], reference["trajectory"])
    window = reference["frames"]
    settings = RgSettings(
        runs=[
            RgRunSettings(label="Protein", selection="protein", calculation_mode="selection"),
            RgRunSettings(
                label="Polymer Oligomers",
                selection="resname SBM EGM",
                calculation_mode="fragments",
                fragment_weighting="equal",
                save_fragment_distribution=True,
                histogram_bins=50,
            ),
        ]
    )
    frames = FrameSelection(start=window["start"], stop=window["stop"], step=window["step"])
    return {item.name: item for item in Rg().compute(universe, frames, settings)}


def test_selection_mode_matches_the_old_per_frame_series(
    reference: dict[str, Any], observables: dict[str, Any]
) -> None:
    """Per-frame protein Rg is unchanged by the port."""
    expected = reference["runs"]["Protein"]["rg_per_frame"]
    np.testing.assert_allclose(
        observables["rg_protein"].values, expected, rtol=0.0, atol=ATOL
    )


def test_fragment_profile_matches_the_old_per_fragment_means(
    reference: dict[str, Any], observables: dict[str, Any]
) -> None:
    """The profile holds the per-fragment mean of the old per-frame matrix."""
    run = reference["runs"]["Polymer Oligomers"]
    profile = observables["rg_polymer_oligomers_fragments"]
    assert len(profile.values) == run["n_fragments"]
    assert profile.index == list(range(run["n_fragments"]))
    np.testing.assert_allclose(
        profile.values, run["fragment_rg_mean"], rtol=0.0, atol=ATOL
    )


def test_fragment_reduction_matches_the_old_equal_weighted_series(
    reference: dict[str, Any], observables: dict[str, Any]
) -> None:
    """The per-frame mean over fragments is unchanged by the port."""
    expected = reference["runs"]["Polymer Oligomers"]["reduced_rg_per_frame_equal"]
    np.testing.assert_allclose(
        observables["rg_polymer_oligomers"].values, expected, rtol=0.0, atol=ATOL
    )


def test_mass_weighted_reduction_matches_the_old_series(reference: dict[str, Any]) -> None:
    """Mass weighting across fragments still reduces the way it used to."""
    universe = _universe(reference["topology"], reference["trajectory"])
    window = reference["frames"]
    settings = RgSettings(
        runs=[
            RgRunSettings(
                label="Polymer Oligomers",
                selection="resname SBM EGM",
                calculation_mode="fragments",
                fragment_weighting="mass",
                save_fragment_distribution=False,
            )
        ]
    )
    frames = FrameSelection(start=window["start"], stop=window["stop"], step=window["step"])
    computed = {item.name: item for item in Rg().compute(universe, frames, settings)}
    np.testing.assert_allclose(
        computed["rg_polymer_oligomers"].values,
        reference["runs"]["Polymer Oligomers"]["reduced_rg_per_frame_mass"],
        rtol=0.0,
        atol=ATOL,
    )


def test_distribution_is_a_profile_over_fixed_bins(observables: dict[str, Any]) -> None:
    """The fragment distribution is a density profile whose index is bin centres.

    The old implementation chose bin edges by pooling every replicate of a
    condition, which a contract plugin cannot do because it sees one replicate.
    Edges now come from ``histogram_range``, so there is nothing to compare
    against the frozen reference; the test checks the shape and the
    normalisation instead.
    """
    distribution = observables["rg_polymer_oligomers_distribution"]
    assert distribution.kind == "profile" and distribution.unit == "1/A"
    assert len(distribution.values) == 50 and len(distribution.index) == 50
    width = distribution.index[1] - distribution.index[0]
    assert sum(distribution.values) * width == pytest.approx(1.0)


def test_fragments_mode_on_a_protein_without_conect_records_raises() -> None:
    """The control topology has no bonds, so fragment mode must refuse to run."""
    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    universe = _universe(CONTROL_TOPOLOGY, CONTROL_TRAJECTORY)
    settings = RgSettings(
        runs=[
            RgRunSettings(label="Protein", selection="protein", calculation_mode="fragments")
        ]
    )
    with pytest.raises(TopologyBondsMissingError):
        Rg().compute(universe, FrameSelection(start=0, stop=2, step=1), settings)


def test_pbc_policy_and_bond_state_reach_the_observable_metadata(
    observables: dict[str, Any],
) -> None:
    """Every observable records how coordinates were treated and where bonds came from."""
    for observable in observables.values():
        assert observable.metadata["pbc_policy"] == "as_loaded"
        assert observable.metadata["topology_has_bonds"] is True
        assert observable.metadata["bond_source"] in {"conect", "guessed"}
