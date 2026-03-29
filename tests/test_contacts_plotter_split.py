"""Tests for contacts plotting functions after extraction to _plotters.py.

Verifies that all 11 private plotting functions are re-exported from
the contacts plugin module (used by plot()), and that 7 shared helpers
exist in the contacts._plotters module (private implementation details).
"""

from __future__ import annotations

import polyzymd.analyses.contacts as contacts_mod
import polyzymd.analyses.contacts._plotters as plotters_mod

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
