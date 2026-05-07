"""Shared constants for analysis and comparison modules.

These constants define default values that are reused across multiple analysis
modules. Centralizing them here ensures consistency and makes them easy to
update.
"""

# Distance cutoff (Angstroms) for polymer-protein contact detection.
DEFAULT_CONTACT_CUTOFF: float = 4.5

# Distance threshold (Angstroms) for catalytic triad / distance pair analyses.
DEFAULT_DISTANCE_THRESHOLD: float = 3.5
