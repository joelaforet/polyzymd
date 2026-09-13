"""Known answers for the contract RMSF plugin.

Each universe here holds four static alpha carbons plus one that hops between
``x = +AMPLITUDE`` and ``x = -AMPLITUDE``. Superposing on the static residues
alone leaves every coordinate untouched, so the fluctuation of the moving
residue is exactly the amplitude and the fluctuation of the others is exactly
zero.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.exceptions import PluginContractError, ReplicateError, SelectionError
from polyzymd.analyses.mda import FrameSelection
from polyzymd.analyses.rmsf import RMSF, RMSFAnalysis, RMSFSettings

AMPLITUDE = 0.5
STATIC = "name CA and resid 1:4"
ALL_CA = "name CA"
N_FRAMES = 6


def make_universe(amplitude: float = AMPLITUDE) -> Any:
    """Five alpha carbons, one of which hops along x between the frames."""
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(5, n_residues=5, atom_resindex=list(range(5)), trajectory=True)
    universe.add_TopologyAttr("names", ["CA"] * 5)
    universe.add_TopologyAttr("resnames", ["ALA"] * 5)
    universe.add_TopologyAttr("resids", [1, 2, 3, 4, 5])
    universe.add_TopologyAttr("masses", [12.0] * 5)
    base = np.asarray(
        [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0], [8.0, 8.0, 8.0]],
        dtype=np.float32,
    )
    frames = []
    for index in range(N_FRAMES):
        positions = base.copy()
        positions[4, 0] += amplitude if index % 2 == 0 else -amplitude
        frames.append(positions)
    universe.load_new(np.stack(frames), format=MemoryReader)
    return universe


def settings(**overrides: Any) -> RMSFSettings:
    """Settings that superpose on the static residues only."""
    fields = {
        "selection": ALL_CA,
        "alignment_selection": STATIC,
        "centroid_selection": STATIC,
    }
    return RMSFSettings(**{**fields, **overrides})


def _frames() -> FrameSelection:
    """The whole synthetic window as a slice."""
    return FrameSelection(start=0, stop=N_FRAMES, step=1, n_frames_total=N_FRAMES)


def profile_of(universe: Any, plugin_settings: RMSFSettings, name: str = "rmsf") -> list[float]:
    """Run one replicate directly and return one observable's values."""
    observables = RMSF().compute(universe, _frames(), plugin_settings)
    return next(observable for observable in observables if observable.name == name).values


def test_only_the_moving_residue_fluctuates() -> None:
    """The profile is the amplitude on residue 5 and zero everywhere else."""
    values = profile_of(make_universe(), settings())

    assert values == pytest.approx([0.0, 0.0, 0.0, 0.0, AMPLITUDE])


def test_the_profile_is_labelled_by_residue_id() -> None:
    """A profile states one residue ID per value, in angstrom."""
    observables = RMSF().compute(make_universe(), _frames(), settings())

    assert len(observables) == 1
    assert observables[0].index == [1.0, 2.0, 3.0, 4.0, 5.0]
    assert observables[0].unit == "A"
    assert observables[0].n_frames == N_FRAMES
    assert observables[0].reduce == "mean_over_index"
    assert observables[0].reduced_kind == "fluctuation"


def test_frame_reference_mode_gives_the_same_fluctuation() -> None:
    """Superposing on static atoms makes the reference choice immaterial here."""
    values = profile_of(make_universe(), settings(reference_mode="frame", reference_frame=1))

    assert values == pytest.approx([0.0, 0.0, 0.0, 0.0, AMPLITUDE])


def test_external_mode_separates_deviation_from_the_reference(tmp_path: Path) -> None:
    """External mode reports the distance from the reference under its own name.

    The reference is frame 0, where the moving residue sits at ``+AMPLITUDE``,
    so its distance from the reference over the window is
    ``sqrt(mean(0, (2 * AMPLITUDE) ** 2))``, which is ``AMPLITUDE * sqrt(2)``,
    while its fluctuation about the trajectory mean stays ``AMPLITUDE``.
    """
    universe = make_universe()
    universe.trajectory[0]
    reference = tmp_path / "reference.pdb"
    universe.atoms.write(str(reference))

    plugin_settings = settings(reference_mode="external", reference_file=str(reference))
    observables = RMSF().compute(make_universe(), _frames(), plugin_settings)

    names = [observable.name for observable in observables]
    assert names == ["rmsf", "rmsd_about_reference_per_residue"]
    assert observables[0].reduced_kind == "fluctuation"
    assert observables[1].reduced_kind is None
    assert observables[0].values == pytest.approx([0.0, 0.0, 0.0, 0.0, AMPLITUDE])
    assert observables[1].values == pytest.approx(
        [0.0, 0.0, 0.0, 0.0, AMPLITUDE * 2**0.5], abs=1e-5
    )


def test_an_empty_selection_is_named() -> None:
    """An empty selection raises rather than reporting zero fluctuation."""
    with pytest.raises(SelectionError, match="matched no atoms"):
        profile_of(make_universe(), settings(selection="name ZZ"))


@pytest.mark.parametrize(
    "overrides",
    [{"reference_mode": "frame"}, {"reference_mode": "external"}],
)
def test_a_reference_mode_without_its_input_is_rejected(overrides: dict[str, Any]) -> None:
    """Frame mode needs a frame number and external mode needs a file."""
    with pytest.raises(ValueError, match="is required when reference_mode"):
        settings(**overrides)


@pytest.mark.parametrize("reference_mode", ["centroid", "average"])
def test_an_explicit_frame_list_is_refused_by_slice_reference_modes(reference_mode: str) -> None:
    """Those modes build their reference from a contiguous slice."""
    frames = FrameSelection(frames=[0, 2, 4], n_frames_total=N_FRAMES)

    with pytest.raises(PluginContractError, match="contiguous trajectory slice"):
        RMSF().compute(make_universe(), frames, settings(reference_mode=reference_mode))


def test_an_explicit_frame_list_is_allowed_by_frame_reference_mode() -> None:
    """Frame mode names its reference, so a non-uniform window is fine."""
    frames = FrameSelection(frames=[0, 2, 4], n_frames_total=N_FRAMES)

    observables = RMSF().compute(
        make_universe(), frames, settings(reference_mode="frame", reference_frame=1)
    )

    assert observables[0].values == pytest.approx([0.0, 0.0, 0.0, 0.0, 0.0])
    assert observables[0].n_frames == 3


def test_an_empty_production_window_is_named() -> None:
    """A window that starts past the last frame is an error, not a zero."""
    frames = FrameSelection(start=N_FRAMES, stop=N_FRAMES + 2, step=1, n_frames_total=N_FRAMES + 2)

    with pytest.raises(ReplicateError, match="holds no frames"):
        RMSF().compute(make_universe(), frames, settings(reference_mode="frame", reference_frame=1))


def test_a_missing_external_reference_is_named(tmp_path: Path) -> None:
    """The file is checked before the trajectory is aligned."""
    missing = tmp_path / "gone.pdb"
    missing.write_text("", encoding="utf-8")
    plugin_settings = settings(reference_mode="external", reference_file=str(missing))
    missing.unlink()

    with pytest.raises(SelectionError, match="does not exist"):
        RMSF().compute(make_universe(), _frames(), plugin_settings)


def test_an_external_reference_with_other_residues_is_named(tmp_path: Path) -> None:
    """A reference holding a different residue set cannot be compared atom by atom."""
    universe = make_universe()
    universe.trajectory[0]
    reference = tmp_path / "partial.pdb"
    universe.select_atoms(STATIC).write(str(reference))

    with pytest.raises(SelectionError, match="same selected atoms"):
        RMSF().compute(
            make_universe(),
            _frames(),
            settings(reference_mode="external", reference_file=str(reference)),
        )


def test_the_external_reference_file_is_part_of_the_identity(tmp_path: Path) -> None:
    """The runner hashes what identity_files names, so the plugin must name it."""
    reference = tmp_path / "reference.pdb"
    reference.write_text("", encoding="utf-8")
    external = settings(reference_mode="external", reference_file=str(reference))

    assert list(RMSF.identity_files(external)) == [reference]
    assert list(RMSF.identity_files(settings())) == []


def test_the_condition_aggregate_carries_the_profile_and_its_mean(
    run_contract_analysis: Any,
) -> None:
    """Three replicates give a per-residue mean and SEM plus the scalar mean."""
    artifact = run_contract_analysis(RMSFAnalysis, settings(), make_universe())

    aggregates = {
        payload["name"]: ObservableAggregate.model_validate(payload)
        for payload in artifact.payload["observables"]
    }
    profile = aggregates["rmsf"]
    scalar = aggregates["rmsf_mean"]

    assert profile.kind == "profile" and profile.n_replicates == 3
    assert profile.profile_mean == pytest.approx([0.0, 0.0, 0.0, 0.0, AMPLITUDE])
    assert profile.profile_sem == pytest.approx([0.0] * 5)
    assert scalar.kind == "fluctuation" and scalar.unit == "A"
    assert scalar.replicate_values == pytest.approx([AMPLITUDE / 5] * 3)
    assert scalar.mean == pytest.approx(AMPLITUDE / 5)
    assert scalar.sem == pytest.approx(0.0)
    assert scalar.n_eff_min is None
