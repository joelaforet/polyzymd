"""Tests for artifact-only contacts plotting helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import polyzymd.analyses.contacts._plotters as plotters_mod
from polyzymd.analyses.mda import ConditionArtifact
from polyzymd.config.comparison import PlotSettings

PLOT_FUNCTIONS = [
    "_plot_contact_fraction_profile",
    "_plot_residence_time_profile",
    "_plot_cf_by_aa_class_bars",
    "_plot_cf_by_partition_bars",
    "_plot_rt_by_aa_class_bars",
    "_plot_rt_by_partition_bars",
]

ARTIFACT_HELPERS = [
    "ContactsProfileData",
    "ContactsConditionPlotData",
    "ContactsPlotData",
    "load_contacts_plot_data",
    "_load_partition_definitions",
]


def test_all_plot_functions_exist_on_plotters() -> None:
    """All stable plot functions should be defined in contacts._plotters."""

    for function_name in PLOT_FUNCTIONS:
        function = getattr(plotters_mod, function_name, None)
        assert function is not None, f"Missing plot function on _plotters: {function_name}"
        assert callable(function), f"{function_name} is not callable"


def test_artifact_plot_helpers_exist() -> None:
    """Artifact-only data containers and loaders should remain available."""

    for helper_name in ARTIFACT_HELPERS:
        helper = getattr(plotters_mod, helper_name, None)
        assert helper is not None, f"Missing artifact helper on _plotters: {helper_name}"


def test_plot_function_count() -> None:
    """There should be exactly six stable contacts plot functions."""

    assert len(PLOT_FUNCTIONS) == 6


def test_artifact_helper_count() -> None:
    """There should be exactly five artifact plotting helpers under test."""

    assert len(ARTIFACT_HELPERS) == 5


def test_rt_by_aa_class_bars_skip_sparse_replicate_dots(tmp_path) -> None:
    """Residence-time AA class bars should not overlay sparse replicate dots."""

    plot_data = _make_plot_data(tmp_path)

    with (
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "rt.png"),
    ):
        paths = plotters_mod._plot_rt_by_aa_class_bars(
            plot_data,
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "rt.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] is None


def test_rt_by_partition_bars_skip_sparse_replicate_dots(tmp_path) -> None:
    """Residence-time partition bars should not overlay sparse replicate dots."""

    plot_data = _make_plot_data(
        tmp_path,
        settings=SimpleNamespace(
            protein_groups={"surface": [1]},
            protein_partitions={"regions": ["surface"]},
        ),
    )

    with (
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "rt_partition.png"),
    ):
        paths = plotters_mod._plot_rt_by_partition_bars(
            plot_data,
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "rt_partition.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] is None


def _make_plot_data(tmp_path: Path, settings: object | None = None):
    """Build a minimal artifact-backed contacts plot dataset."""

    profile = plotters_mod.ContactsProfileData(
        replicate_ids=np.asarray([1, 2], dtype=np.int64),
        residue_ids=np.asarray([1], dtype=np.int64),
        residue_names=np.asarray(["ALA"], dtype=str),
        residue_groups=np.asarray(["nonpolar"], dtype=str),
        contact_fraction_by_replicate=np.asarray([[0.4], [0.6]], dtype=np.float64),
        contact_fraction_mean=np.asarray([0.5], dtype=np.float64),
        contact_fraction_sem=np.asarray([0.1], dtype=np.float64),
        polymer_types=np.asarray(["PEG"], dtype=str),
        contact_fraction_by_polymer_type=np.asarray([[[0.4], [0.6]]], dtype=np.float64),
        residence_time_mean_ns=np.asarray([[5.0]], dtype=np.float64),
        residence_time_sem_ns=np.asarray([[0.5]], dtype=np.float64),
        residence_time_event_counts=np.asarray([[2]], dtype=np.int64),
    )
    artifact = ConditionArtifact(
        analysis_name="contacts",
        condition_label="Control",
        replicates=[1, 2],
        payload={},
        metadata={"compute_residence_times": True},
    )
    condition = plotters_mod.ContactsConditionPlotData(
        label="Control",
        aggregated_dir=tmp_path / "aggregated",
        artifact=artifact,
        profile=profile,
    )
    return plotters_mod.ContactsPlotData(
        conditions={"Control": condition},
        labels=("Control",),
        settings=settings or SimpleNamespace(protein_groups=None, protein_partitions=None),
    )
