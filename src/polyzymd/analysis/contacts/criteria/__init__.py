"""Contact criteria using the Strategy pattern.

This module provides different methods for determining whether
two molecular groups are "in contact". Users can define custom
criteria by subclassing ContactCriteria.

Built-in criteria:
- AnyAtomWithinCutoff: Contact if any atom pair within cutoff
- AnyAtomToCOM: Contact if any atom within cutoff of target COM
- COMToCOM: Contact if COMs within cutoff
- MinimumDistance: Returns minimum distance (for continuous analysis)

Examples
--------
>>> # Default: any atom within 4 Angstroms
>>> criteria = AnyAtomWithinCutoff(cutoff=4.0)
>>>
>>> # Custom: polymer atom to protein residue COM
>>> criteria = AnyAtomToCOM(cutoff=4.0)
>>>
>>> # Define your own criteria
>>> class HydrogenBondCriteria(ContactCriteria):
...     def is_contact(self, query_atoms, target_atoms, box=None):
...         # Custom hydrogen bond detection logic
...         ...
"""

from polyzymd.analysis.contacts.criteria.base import ContactCriteria
from polyzymd.analysis.contacts.criteria.distance import (
    AnyAtomWithinCutoff,
    AnyAtomToCOM,
    COMToCOM,
    MinimumDistance,
)

__all__ = [
    "ContactCriteria",
    "AnyAtomWithinCutoff",
    "AnyAtomToCOM",
    "COMToCOM",
    "MinimumDistance",
]
