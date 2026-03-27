"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._contacts_results``.

All public names are re-exported so that existing ``from polyzymd.analysis.contacts.results
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._contacts_results import (  # noqa: F401
    ContactEvent,
    ContactResult,
    PolymerSegmentContacts,
    ResidueContactData,
    compress_contact_array,
    decompress_events,
)

__all__ = [
    "ContactEvent",
    "ContactResult",
    "PolymerSegmentContacts",
    "ResidueContactData",
    "compress_contact_array",
    "decompress_events",
]
