"""Generic modifier aliases for conjugation construction workflows."""

from __future__ import annotations

from polyzymd.builders.conjugation.pdb_assembly import PlacedPolymerFragment
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment

GeneratedModifierFragment = GeneratedPolymerFragment
PlacedModifierFragment = PlacedPolymerFragment

__all__ = ["GeneratedModifierFragment", "PlacedModifierFragment"]
