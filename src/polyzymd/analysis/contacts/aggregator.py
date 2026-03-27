"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._contacts_aggregator``.

All public names are re-exported so that existing ``from polyzymd.analysis.contacts.aggregator
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._contacts_aggregator import (  # noqa: F401
    AggregatedContactResult,
    AggregatedResidueStats,
    aggregate_contact_results,
    compute_mad,
    compute_sem,
    load_and_aggregate,
)

__all__ = [
    "AggregatedContactResult",
    "AggregatedResidueStats",
    "aggregate_contact_results",
    "compute_mad",
    "compute_sem",
    "load_and_aggregate",
]
