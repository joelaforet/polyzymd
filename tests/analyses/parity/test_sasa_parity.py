"""Parity of the contract SASA plugin with the implementation it replaced.

``sasa_reference.json`` holds per-frame total areas and per-residue mean areas
frozen from the pre-port ``polyzymd.analyses.sasa._mda`` kernel, on frames 500
to 539 of two real LipA trajectories, before that code was deleted. The file
records the topology and trajectory paths, the frame window, the settings and
the commit it was frozen at.

The port reuses the same MDTraj call and the same nm^2 to A^2 conversion in the
same order, so the values must agree to floating-point noise. The tolerance is
an absolute 1e-6 A^2 on totals of roughly 9000 A^2, which is far tighter than
any physically meaningful difference and loose enough to survive a change in
summation order inside NumPy.

The tests skip when the trajectories are not on this machine, so a checkout
without the data still passes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.sasa import SASA, SASARun, SASASettings

REFERENCE = Path(__file__).with_name("sasa_reference.json")
TOLERANCE_A2 = 1e-6


def _reference() -> dict[str, Any]:
    """Frozen values and the provenance block that describes them."""
    return json.loads(REFERENCE.read_text())


def _universe(entry: dict[str, Any]) -> Any:
    """Universe for one frozen replicate, or a skip when the files are absent.

    The topology md5 is checked before anything is measured, so a reference
    that no longer describes the file on disk fails loudly instead of failing
    as a wrong number.
    """
    import MDAnalysis as mda

    topology, trajectory = Path(entry["topology"]), Path(entry["trajectory"])
    if not (topology.exists() and trajectory.exists()):
        pytest.skip(f"real trajectory not on this machine: {trajectory}")
    digest = hashlib.md5(topology.read_bytes()).hexdigest()
    assert digest == entry["topology_md5"], f"{topology} is not the frozen topology"
    universe = mda.Universe(str(topology), str(trajectory))
    assert universe.atoms.n_atoms == entry["n_atoms"]
    return universe


@pytest.mark.parametrize("dataset", ["control", "polymer_50_50"])
def test_sasa_matches_the_frozen_reference(dataset: str) -> None:
    """Every total and every per-residue relative area reproduces the old kernel."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    reference = _reference()
    provenance, entry = reference["provenance"], reference["datasets"][dataset]
    universe = _universe(entry)
    window = provenance["frames"]
    settings = SASASettings(
        runs=[
            SASARun(
                label=label,
                target_selection=run["target_selection"],
                context_selection=run["context_selection"],
            )
            for label, run in entry["runs"].items()
        ],
        **provenance["settings"],
    )

    observables = {
        observable.name: observable
        for observable in SASA().compute(
            universe,
            FrameSelection(start=window["start"], stop=window["stop"], step=window["step"]),
            settings,
        )
    }

    for label, run in entry["runs"].items():
        total = observables[f"sasa_{label}"]
        profile = observables[f"relative_sasa_{label}"]
        assert total.unit == "A^2"
        assert len(total.values) == window["stop"] - window["start"]
        np.testing.assert_allclose(total.values, run["total_sasa_a2"], rtol=0, atol=TOLERANCE_A2)
        labels = profile.metadata["residue_labels"]
        assert [label.split(":")[1] for label in labels] == [
            str(resid) for resid in run["residue_resids"]
        ]
        assert [label.split(":")[2] for label in labels] == run["residue_resnames"]
        assert profile.index_label == "residue index"
        np.testing.assert_allclose(
            profile.values, run["residue_mean_relative_sasa"], rtol=0, atol=TOLERANCE_A2
        )


def test_reference_records_its_provenance() -> None:
    """The frozen file names the code, the window and the settings it came from."""
    provenance = _reference()["provenance"]

    assert "pre-port" in provenance["source"]
    assert provenance["polyzymd_commit"]
    assert provenance["settings"] == {
        "probe_radius_nm": 0.14,
        "n_sphere_points": 960,
        "chunk_size": 50,
    }
    for entry in _reference()["datasets"].values():
        assert len(entry["topology_md5"]) == 32
