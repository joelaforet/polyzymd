"""Tests for contacts plotting functions after extraction to _plotters.py.

Verifies that all 11 private plotting functions are re-exported from
the contacts plugin module (used by plot()), and that 7 shared helpers
exist in the contacts._plotters module (private implementation details).
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import polyzymd.analyses.contacts as contacts_mod
import polyzymd.analyses.contacts._plotters as plotters_mod
from polyzymd.config.comparison import PlotSettings

# All 11 plot functions — re-exported from __init__.py (used in plot())
PLOT_FUNCTIONS = [
    "_plot_contact_fraction_profile",
    "_plot_residence_time_profile",
    "_plot_cf_by_aa_class_bars",
    "_plot_cf_by_partition_bars",
    "_plot_rt_by_aa_class_bars",
    "_plot_rt_by_partition_bars",
    "_plot_user_partition_bars",
    "_plot_system_coverage_bars",
    "_plot_system_coverage_heatmap",
    "_plot_binding_preference_bars",
    "_plot_binding_preference_heatmap",
]

# All 7 shared helpers — private to _plotters.py (data loaders + utilities)
SHARED_HELPERS = [
    "_get_polymer_types_and_aa_classes",
    "_get_enrichment_value",
    "_get_enrichment_with_sem",
    "_load_binding_preference_results",
    "_load_system_coverage_results",
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
    """There should be exactly 11 plot functions."""
    assert len(PLOT_FUNCTIONS) == 11


def test_shared_helper_count() -> None:
    """There should be exactly 7 shared helpers."""
    assert len(SHARED_HELPERS) == 7


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


def test_system_coverage_bars_pass_replicate_enrichments(tmp_path) -> None:
    """System coverage bars should pass aligned per-replicate enrichments."""
    aromatic = MagicMock(
        mean_coverage_enrichment=0.5,
        sem_coverage_enrichment=0.1,
        per_replicate_enrichments=[0.4, 0.6],
    )
    polar = MagicMock(
        mean_coverage_enrichment=-0.2,
        sem_coverage_enrichment=0.05,
        per_replicate_enrichments=[-0.3, -0.1],
    )
    coverage = MagicMock()
    coverage.aa_class_names.return_value = ["aromatic", "polar"]
    coverage.aa_class_coverage.get_entry.side_effect = {
        "aromatic": aromatic,
        "polar": polar,
    }.get

    with (
        patch.object(
            plotters_mod, "_load_system_coverage_results", return_value={"CondA": coverage}
        ),
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "coverage.png"),
    ):
        paths = plotters_mod._plot_system_coverage_bars(
            {},
            ["CondA"],
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "coverage.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] == [[[0.4, 0.6], [-0.3, -0.1]]]


def test_user_partition_bars_pass_replicate_enrichments(tmp_path) -> None:
    """User partition bars should align replicate enrichments by element."""
    lid = MagicMock(
        mean_coverage_enrichment=1.2,
        sem_coverage_enrichment=0.2,
        per_replicate_enrichments=[1.0, 1.4],
    )
    core = MagicMock(
        mean_coverage_enrichment=-0.4,
        sem_coverage_enrichment=0.1,
        per_replicate_enrichments=[-0.5, -0.3],
    )
    partition = MagicMock()
    partition.element_names.return_value = ["lid", "core"]
    partition.get_entry.side_effect = {"lid": lid, "core": core}.get
    coverage = MagicMock()
    coverage.user_defined_partitions = {"regions": partition}

    with (
        patch.object(
            plotters_mod, "_load_system_coverage_results", return_value={"CondA": coverage}
        ),
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "regions.png"),
    ):
        paths = plotters_mod._plot_user_partition_bars(
            {},
            ["CondA"],
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "regions.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] == [[[1.0, 1.4], [-0.5, -0.3]]]


def test_binding_preference_bars_pass_replicate_enrichments(tmp_path) -> None:
    """Binding preference bars should align replicate enrichments by group."""
    aromatic = MagicMock(
        mean_enrichment=0.8,
        sem_enrichment=0.2,
        per_replicate_enrichments=[0.6, 1.0],
    )
    polar = MagicMock(
        mean_enrichment=-0.1,
        sem_enrichment=0.05,
        per_replicate_enrichments=[-0.2, 0.0],
    )
    result = MagicMock()
    result.binding_preference = None
    result.get_entry.side_effect = lambda polymer_type, group: {
        ("PEG", "aromatic"): aromatic,
        ("PEG", "polar"): polar,
    }.get((polymer_type, group))

    with (
        patch.object(
            plotters_mod, "_load_binding_preference_results", return_value={"CondA": result}
        ),
        patch.object(
            plotters_mod,
            "_get_polymer_types_and_aa_classes",
            return_value=(["PEG"], ["aromatic", "polar"]),
        ),
        patch.object(plotters_mod, "grouped_bars") as mock_grouped_bars,
        patch.object(plotters_mod, "apply_legend"),
        patch.object(plotters_mod, "save_figure", return_value=tmp_path / "binding.png"),
    ):
        paths = plotters_mod._plot_binding_preference_bars(
            {},
            ["CondA"],
            tmp_path,
            PlotSettings(),
        )

    assert paths == [tmp_path / "binding.png"]
    assert mock_grouped_bars.call_args.kwargs["replicate_values"] == [[[0.6, 1.0], [-0.2, 0.0]]]
