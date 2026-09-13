"""Known-answer tests for the rg contract plugin.

The synthetic universes here have radii of gyration that can be worked out by
hand: four unit-mass atoms on a cross of scale s have a radius of gyration of
exactly s about their centre of mass.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.exceptions import (
    ReplicateError,
    SelectionError,
    TopologyBondsMissingError,
)
from polyzymd.analyses.rg import Rg, RgAnalysis, RgRunSettings, RgSettings
from tests.analyses.conftest import CROSS, make_synthetic_universe

FRAMES = (0, 3, 1)


def _window() -> Any:
    """Frame selection covering the whole synthetic trajectory."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    start, stop, step = FRAMES
    return FrameSelection(start=start, stop=stop, step=step)


def _bonded_universe(
    scales: tuple[float, ...] = (1.0, 2.0), masses: tuple[float, ...] = (1.0, 3.0)
):
    """Universe of one bonded cross per scale, so fragment mode has fragments."""
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    n_fragments = len(scales)
    universe = mda.Universe.empty(
        4 * n_fragments,
        n_residues=n_fragments,
        atom_resindex=np.repeat(np.arange(n_fragments), 4),
        trajectory=True,
    )
    universe.add_TopologyAttr("masses", np.repeat(masses, 4).tolist())
    universe.add_TopologyAttr("resnames", ["POL"] * n_fragments)
    universe.add_TopologyAttr(
        "bonds", [(4 * i + 0, 4 * i + j) for i in range(n_fragments) for j in (1, 2, 3)]
    )
    positions = np.concatenate(
        [
            np.asarray(CROSS, dtype=np.float64) * scale + (30.0 * i, 0.0, 0.0)
            for i, scale in enumerate(scales)
        ]
    )
    universe.load_new(np.stack([positions] * 3).astype(np.float32), format=MemoryReader)
    return universe


def _compute(universe: Any, *runs: RgRunSettings) -> dict[str, Any]:
    """Run the plugin over the whole trajectory and key the result by name."""
    return {
        observable.name: observable
        for observable in Rg().compute(universe, _window(), RgSettings(runs=list(runs)))
    }


def test_selection_mode_measures_the_whole_group() -> None:
    """A cross of scale 2 has a radius of gyration of exactly 2."""
    observables = _compute(
        make_synthetic_universe(scale=2.0, n_frames=3),
        RgRunSettings(label="Protein", selection="all"),
    )

    observable = observables["rg_protein"]
    assert observable.kind == "mean_of_timeseries" and observable.unit == "A"
    assert observable.values == pytest.approx([2.0, 2.0, 2.0])
    assert observable.metadata["pbc_policy"] == "as_loaded"


def test_a_label_with_spaces_becomes_one_slug() -> None:
    """The observable name is the label slugged, so YAML labels can be prose."""
    observables = _compute(
        make_synthetic_universe(), RgRunSettings(label="Polymer Oligomers", selection="all")
    )

    assert set(observables) == {"rg_polymer_oligomers"}


def test_fragments_mode_reports_the_profile_and_its_mean() -> None:
    """Two bonded crosses of scale 1 and 2 give a profile of [1, 2] and a mean of 1.5."""
    observables = _compute(
        _bonded_universe(),
        RgRunSettings(
            label="polymer",
            selection="all",
            calculation_mode="fragments",
            histogram_range=(0.0, 5.0),
        ),
    )

    profile = observables["rg_polymer_fragments"]
    assert profile.kind == "profile" and profile.index == [0.0, 1.0]
    assert profile.values == pytest.approx([1.0, 2.0])
    assert observables["rg_polymer"].values == pytest.approx([1.5, 1.5, 1.5])


def test_mass_weighting_leans_on_the_heavier_fragment() -> None:
    """Fragment masses of 4 and 12 pull the mean of 1 and 2 to 1.75."""
    observables = _compute(
        _bonded_universe(),
        RgRunSettings(
            label="polymer",
            selection="all",
            calculation_mode="fragments",
            fragment_weighting="mass",
            save_fragment_distribution=False,
        ),
    )

    assert observables["rg_polymer"].values == pytest.approx([1.75, 1.75, 1.75])


def test_the_distribution_is_a_density_over_bin_centres() -> None:
    """The distribution profile integrates to one over the histogram range."""
    observables = _compute(
        _bonded_universe(),
        RgRunSettings(
            label="polymer",
            selection="all",
            calculation_mode="fragments",
            histogram_bins=10,
            histogram_range=(0.0, 5.0),
        ),
    )

    distribution = observables["rg_polymer_distribution"]
    assert distribution.index == pytest.approx([0.25 + 0.5 * i for i in range(10)])
    assert sum(distribution.values) * 0.5 == pytest.approx(1.0)


def test_a_fragment_outside_the_histogram_range_raises() -> None:
    """Silently dropping a fragment from the density would misreport the shape."""
    run = RgRunSettings(
        label="polymer",
        selection="all",
        calculation_mode="fragments",
        histogram_range=(0.0, 1.5),
    )

    with pytest.raises(ReplicateError, match="outside histogram_range"):
        _compute(_bonded_universe(), run)


def test_fragment_mode_without_bonds_raises() -> None:
    """A topology with no bonds makes fragment mode meaningless, so it fails."""
    run = RgRunSettings(
        label="polymer",
        selection="all",
        calculation_mode="fragments",
        save_fragment_distribution=False,
    )

    with pytest.raises(TopologyBondsMissingError):
        _compute(make_synthetic_universe(), run)


def test_the_single_fragment_fallback_is_still_available() -> None:
    """Opting in measures the whole selection as one fragment, as before 1.3."""
    observables = _compute(
        make_synthetic_universe(scale=2.0, n_frames=3),
        RgRunSettings(
            label="polymer",
            selection="all",
            calculation_mode="fragments",
            allow_single_fragment_fallback=True,
            save_fragment_distribution=False,
        ),
    )

    assert observables["rg_polymer_fragments"].values == pytest.approx([2.0])
    assert observables["rg_polymer"].values == pytest.approx([2.0, 2.0, 2.0])


def test_an_empty_selection_raises_instead_of_reporting_zero() -> None:
    """An empty selection is a configuration error, not a radius of gyration of 0."""
    run = RgRunSettings(label="protein", selection="index 99")

    with pytest.raises(SelectionError, match="matched no atoms"):
        _compute(make_synthetic_universe(), run)


def test_fragment_weighting_is_rejected_on_a_selection_run() -> None:
    """A setting that does nothing in the chosen mode is a mistake worth naming."""
    with pytest.raises(ValueError, match="fragment_weighting"):
        RgRunSettings(label="protein", selection="all", fragment_weighting="mass")


def test_two_labels_that_slug_the_same_are_rejected() -> None:
    """Two runs may not name the same observable."""
    with pytest.raises(ValueError, match="unique after slugging"):
        RgSettings(
            runs=[
                RgRunSettings(label="Protein A", selection="all"),
                RgRunSettings(label="protein a", selection="all"),
            ]
        )


CAMPAIGN_RUNS = [
    {"label": "Protein", "selection": "protein", "calculation_mode": "selection"},
    {
        "label": "Polymer Oligomers",
        "selection": "resname SBM EGM",
        "calculation_mode": "fragments",
        "fragment_weighting": "equal",
        "save_fragment_distribution": True,
        "histogram_bins": 50,
    },
]


def test_the_campaign_settings_need_only_the_bin_range_added() -> None:
    """Every other key of the LipA comparison file still parses as it stands."""
    runs = [dict(run) for run in CAMPAIGN_RUNS]
    runs[1]["histogram_range"] = [6.0, 10.0]

    settings = RgSettings.model_validate({"runs": runs})

    assert [run.label for run in settings.runs] == ["Protein", "Polymer Oligomers"]
    assert settings.runs[1].histogram_range == (6.0, 10.0)


def test_a_distribution_without_a_range_is_rejected_by_name() -> None:
    """The bins cannot be guessed, so the error says which setting is missing and why."""
    with pytest.raises(ValueError, match="histogram_range") as excinfo:
        RgSettings.model_validate({"runs": [dict(run) for run in CAMPAIGN_RUNS]})

    message = str(excinfo.value)
    assert "save_fragment_distribution" in message
    assert "every replicate" in message


def test_the_metadata_reaches_the_written_replicate_artifact(
    tmp_path: Path, run_contract_analysis: Callable[..., Any]
) -> None:
    """The framework carries per-observable provenance onto disk, not just in memory."""
    from polyzymd.analyses.mda import ArtifactStore

    run_contract_analysis(
        RgAnalysis,
        RgSettings(runs=[RgRunSettings(label="Protein", selection="all")]),
        make_synthetic_universe(),
        root=tmp_path,
    )

    artifact = ArtifactStore(tmp_path / "analysis" / "A" / "rg" / "run_1").read_replicate_result()
    metadata = artifact.payload["observables"][0]["metadata"]
    assert metadata["pbc_policy"] == "as_loaded"
    assert metadata["topology_has_bonds"] is False
    assert metadata["bond_source"] == "none"


def test_the_lifecycle_aggregates_three_replicates(
    run_contract_analysis: Callable[..., Any],
) -> None:
    """Three replicates of scale 1.5, 2.0 and 2.5 average to 2.0."""
    artifact = run_contract_analysis(
        RgAnalysis,
        RgSettings(runs=[RgRunSettings(label="Protein", selection="all")]),
        lambda replicate: make_synthetic_universe(scale=1.0 + 0.5 * replicate),
    )

    aggregate = ObservableAggregate.model_validate(artifact.payload["observables"][0])
    assert aggregate.name == "rg_protein" and aggregate.n_replicates == 3
    assert aggregate.replicate_values == pytest.approx([1.5, 2.0, 2.5], abs=1e-5)
    assert aggregate.mean == pytest.approx(2.0, abs=1e-5)
