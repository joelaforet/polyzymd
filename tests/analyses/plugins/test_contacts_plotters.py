"""Tests for contacts plotting functions after extraction to _plotters.py.

Verifies that stable private plotting functions are re-exported from the
contacts plugin module and that shared contact plotting helpers remain
available in the contacts._plotters module.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import polyzymd.analyses.contacts as contacts_mod
import polyzymd.analyses.contacts._plotters as plotters_mod
from polyzymd.config.comparison import PlotSettings

# Stable plot functions — re-exported from __init__.py (used in plot())
PLOT_FUNCTIONS = [
    "_plot_contact_fraction_profile",
    "_plot_residence_time_profile",
    "_plot_cf_by_aa_class_bars",
    "_plot_cf_by_partition_bars",
    "_plot_rt_by_aa_class_bars",
    "_plot_rt_by_partition_bars",
]

# Shared helpers — private to _plotters.py (data loaders + utilities)
SHARED_HELPERS = [
    "_load_aggregated_contact_results",
    "_load_partition_definitions",
]


def test_all_plot_functions_exist_on_module() -> None:
    """All 11 plot functions should be re-exported from contacts __init__.py."""
    for fn_name in PLOT_FUNCTIONS:
        fn = getattr(contacts_mod, fn_name, None)
        assert fn is not None, f"Missing plot function on contacts module: {fn_name}"
        assert callable(fn), f"{fn_name} is not callable"


def test_all_plot_functions_exist_on_plotters() -> None:
    """All 11 plot functions should be defined in contacts._plotters."""
    for fn_name in PLOT_FUNCTIONS:
        fn = getattr(plotters_mod, fn_name, None)
        assert fn is not None, f"Missing plot function on _plotters: {fn_name}"
        assert callable(fn), f"{fn_name} is not callable"


def test_all_shared_helpers_exist() -> None:
    """All 7 shared helpers should be callable in contacts._plotters."""
    for fn_name in SHARED_HELPERS:
        fn = getattr(plotters_mod, fn_name, None)
        assert fn is not None, f"Missing shared helper on _plotters: {fn_name}"
        assert callable(fn), f"{fn_name} is not callable"


def test_plot_function_count() -> None:
    """There should be exactly 6 stable plot functions."""
    assert len(PLOT_FUNCTIONS) == 6


def test_shared_helper_count() -> None:
    """There should be exactly 2 shared helpers."""
    assert len(SHARED_HELPERS) == 2


def test_rt_by_aa_class_bars_skip_sparse_replicate_dots(tmp_path) -> None:
    """Residence-time AA class bars should not overlay sparse replicate dots."""
    result = MagicMock()
    result.residue_stats = [
        MagicMock(
            protein_resid=1,
            protein_group="nonpolar",
            residence_time_by_polymer_type={"PEG": (5.0, 0.5)},
        )
    ]
    result.polymer_types.return_value = {"PEG"}
    result.group_residence_time.return_value = {"nonpolar": (5.0, 0.5)}

    with (
        patch.object(
            plotters_mod, "_load_aggregated_contact_results", return_value={"Control": result}
        ),
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "rt.png"),
    ):
        paths = plotters_mod._plot_rt_by_aa_class_bars(
            {},
            ["Control"],
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "rt.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] is None
    result.group_residence_time_per_replicate.assert_not_called()


def test_rt_by_partition_bars_skip_sparse_replicate_dots(tmp_path) -> None:
    """Residence-time partition bars should not overlay sparse replicate dots."""
    result = MagicMock()
    result.residue_stats = [
        MagicMock(
            protein_resid=1,
            protein_group="nonpolar",
            residence_time_by_polymer_type={"PEG": (5.0, 0.5)},
        )
    ]
    result.polymer_types.return_value = {"PEG"}
    result.subset_residence_time.return_value = (5.0, 0.5)

    with (
        patch.object(
            plotters_mod, "_load_aggregated_contact_results", return_value={"Control": result}
        ),
        patch.object(
            plotters_mod,
            "_load_partition_definitions",
            return_value=({"surface": [1]}, {"regions": ["surface"]}),
        ),
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "rt_partition.png"),
    ):
        paths = plotters_mod._plot_rt_by_partition_bars(
            {},
            ["Control"],
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "rt_partition.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] is None
    result.subset_residence_time_per_replicate.assert_not_called()
