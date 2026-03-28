"""Tests for inlined contacts plotting functions.

Verifies that all 11 private plotting functions and 7 shared helpers
exist as module-level callables in the contacts plugin module.
"""

from __future__ import annotations

import polyzymd.analyses.contacts as contacts_mod

# All 11 inlined plot functions
INLINED_PLOT_FUNCTIONS = [
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

# All 7 inlined shared helpers
INLINED_SHARED_HELPERS = [
    "_get_polymer_types_and_aa_classes",
    "_get_enrichment_value",
    "_get_enrichment_with_sem",
    "_load_binding_preference_results",
    "_load_system_coverage_results",
    "_load_aggregated_contact_results",
    "_load_partition_definitions",
]


def test_all_plot_functions_exist() -> None:
    """All 11 inlined plot functions should be callable attributes of the module."""
    for fn_name in INLINED_PLOT_FUNCTIONS:
        fn = getattr(contacts_mod, fn_name, None)
        assert fn is not None, f"Missing plot function: {fn_name}"
        assert callable(fn), f"{fn_name} is not callable"


def test_all_shared_helpers_exist() -> None:
    """All 7 inlined shared helpers should be callable attributes of the module."""
    for fn_name in INLINED_SHARED_HELPERS:
        fn = getattr(contacts_mod, fn_name, None)
        assert fn is not None, f"Missing shared helper: {fn_name}"
        assert callable(fn), f"{fn_name} is not callable"


def test_plot_function_count() -> None:
    """There should be exactly 11 inlined plot functions."""
    assert len(INLINED_PLOT_FUNCTIONS) == 11


def test_shared_helper_count() -> None:
    """There should be exactly 7 inlined shared helpers."""
    assert len(INLINED_SHARED_HELPERS) == 7
