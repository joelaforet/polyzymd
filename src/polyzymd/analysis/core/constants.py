"""Shared constants for analysis and comparison modules.

These constants define default values that are reused across multiple modules
(analysis configs, comparison settings, surface exposure filters, helpers).
Centralizing them here ensures consistency and makes them easy to update.
"""

# Relative SASA threshold for classifying a residue as surface-exposed.
# Residues with SASA/maxSASA > threshold are considered accessible.
# Used in: ContactsConfig, ContactsAnalysisSettings, ExposureAnalysisSettings,
#          BindingFreeEnergyAnalysisSettings, SurfaceExposureResult,
#          SurfaceExposureFilter, compute_condition_binding_preference().
DEFAULT_SURFACE_EXPOSURE_THRESHOLD: float = 0.2

# Distance cutoff (Angstroms) for polymer-protein contact detection.
# Used in: ContactsConfig, ContactsAnalysisSettings.
DEFAULT_CONTACT_CUTOFF: float = 4.5

# Distance threshold (Angstroms) for catalytic triad / distance pair analyses.
# Used in: CatalyticTriadConfig, DistancesAnalysisSettings.
DEFAULT_DISTANCE_THRESHOLD: float = 3.5
