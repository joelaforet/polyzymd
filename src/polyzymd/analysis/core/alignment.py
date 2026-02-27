"""Trajectory alignment utilities.

This module provides reusable trajectory alignment infrastructure for
all analysis modules. The alignment is performed in-memory using
MDAnalysis and removes rotational drift and center-of-mass motion
that can add noise to distance and flexibility measurements.

Reference Modes
---------------
- **centroid**: Align to most populated conformational cluster (K-Means).
  Best for measuring fluctuations around the equilibrium state.
- **average**: Align to mathematical average structure.
  Best for pure thermal fluctuation analysis (note: average may have
  unphysical geometry).
- **frame**: Align to a specific user-specified frame number.
  Best for comparing to a known functional conformation.
- **external**: Align to an external PDB structure (e.g., crystal structure).
  Best for measuring deviations from a catalytically competent reference
  geometry that is independent of the simulation conditions.

Usage
-----
>>> from polyzymd.analysis.core.alignment import AlignmentConfig, align_trajectory
>>>
>>> # Create configuration
>>> config = AlignmentConfig(
...     enabled=True,
...     reference_mode="centroid",
...     selection="protein and name CA",
... )
>>>
>>> # Align trajectory in-place
>>> ref_frame = align_trajectory(universe, config, start_frame=100)
>>> # Universe is now aligned; ref_frame is the reference frame index

Notes
-----
Alignment modifies the Universe in-place. If you need the original
coordinates, load a fresh Universe after alignment.

The alignment selection should typically be stable atoms that define
the reference frame. Common choices:
- "protein and name CA": Backbone alpha-carbons (default)
- "protein and backbone": Full backbone (N, CA, C)
- "protein": All protein atoms (includes flexible side chains)

References
----------
- MDAnalysis alignment: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from pydantic import BaseModel, Field, model_validator

from polyzymd.analysis.core.loader import _require_mdanalysis

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

LOGGER = logging.getLogger(__name__)

# Type alias for reference modes (re-export from centroid for consistency)
ReferenceMode = Literal["centroid", "average", "frame", "external"]


class AlignmentConfig(BaseModel):
    """Configuration for trajectory alignment.

    This model defines how a trajectory should be aligned before analysis.
    By default, alignment is ENABLED to remove rotational drift and COM motion
    that can add noise to distance measurements.

    Attributes
    ----------
    enabled : bool
        Whether to align the trajectory before analysis. Default is True.
        When enabled, the user is notified via INFO-level logging.
    reference_mode : ReferenceMode
        How to select the reference structure:
        - "centroid": Most populated cluster center (K-Means)
        - "average": Mathematical average structure
        - "frame": Specific frame number
    reference_frame : int | None
        Frame number (1-indexed) when reference_mode="frame".
        Required when mode is "frame", ignored otherwise.
    selection : str
        MDAnalysis selection string for superposition. Atoms matching
        this selection are used to compute the optimal rotation/translation.
        Default: "protein and name CA" (backbone alpha-carbons).
    centroid_selection : str
        Selection for K-Means clustering when mode="centroid".
        Default: "protein" (all protein atoms for better discrimination).

    Examples
    --------
    >>> # Default: align to centroid using CA atoms
    >>> config = AlignmentConfig()

    >>> # Align to a specific frame
    >>> config = AlignmentConfig(
    ...     reference_mode="frame",
    ...     reference_frame=500,
    ... )

    >>> # Disable alignment (not recommended for most analyses)
    >>> config = AlignmentConfig(enabled=False)
    """

    enabled: bool = Field(
        default=True,
        description="Whether to align trajectory before analysis",
    )
    reference_mode: ReferenceMode = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, or frame",
    )
    reference_frame: int | None = Field(
        default=None,
        description="Frame number (1-indexed) when reference_mode='frame'",
    )
    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for superposition",
    )
    centroid_selection: str = Field(
        default="protein",
        description="Selection for K-Means clustering (centroid mode only)",
    )
    reference_file: Path | None = Field(
        default=None,
        description=(
            "Path to an external PDB file to use as the alignment reference "
            "(required when reference_mode='external'). The external PDB "
            "must contain protein atoms that match the simulation topology."
        ),
    )

    @model_validator(mode="after")
    def validate_reference_params(self) -> "AlignmentConfig":
        """Validate reference_frame and reference_file for their modes."""
        if not self.enabled:
            return self
        if self.reference_mode == "frame" and self.reference_frame is None:
            raise ValueError(
                "reference_frame is required when reference_mode='frame'. "
                "Provide a 1-indexed frame number."
            )
        if self.reference_mode == "external":
            if self.reference_file is None:
                raise ValueError(
                    "reference_file is required when reference_mode='external'. "
                    "Provide a path to the external PDB reference structure."
                )
            ref_path = Path(self.reference_file)
            if not ref_path.exists():
                raise ValueError(
                    f"reference_file does not exist: {ref_path}. "
                    "Provide a valid path to the external PDB reference structure."
                )
        return self

    def to_dict(self) -> dict:
        """Convert to dictionary for cache key hashing."""
        d = {
            "enabled": self.enabled,
            "reference_mode": self.reference_mode,
            "reference_frame": self.reference_frame,
            "selection": self.selection,
            "centroid_selection": self.centroid_selection,
        }
        if self.reference_file is not None:
            d["reference_file"] = str(self.reference_file)
        return d


def align_trajectory(
    universe: "Universe",
    config: AlignmentConfig,
    start_frame: int = 0,
    stop_frame: int | None = None,
) -> int | None:
    """Align trajectory in-memory to a reference structure.

    This function performs trajectory alignment using MDAnalysis, removing
    rotational drift and center-of-mass motion. The alignment is performed
    in-place, modifying the Universe's coordinates.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe to align. Will be modified in-place.
    config : AlignmentConfig
        Alignment configuration specifying reference mode and selections.
    start_frame : int, optional
        First frame to consider (0-indexed). Default is 0.
        Frames before this are typically equilibration and should be excluded.
    stop_frame : int, optional
        Last frame (exclusive). If None, uses all frames after start_frame.

    Returns
    -------
    int | None
        Reference frame index (0-indexed), or None if using average mode.
        This can be stored in results for reproducibility.

    Raises
    ------
    ImportError
        If MDAnalysis is not available.
    ValueError
        If reference_mode is "frame" but the frame number is out of range.

    Notes
    -----
    When alignment is performed, an INFO-level log message is emitted to
    notify the user. This ensures users are aware that their trajectory
    coordinates have been modified.

    Examples
    --------
    >>> from polyzymd.analysis.core.alignment import AlignmentConfig, align_trajectory
    >>>
    >>> config = AlignmentConfig(reference_mode="centroid")
    >>> ref_frame = align_trajectory(universe, config, start_frame=100)
    >>> print(f"Aligned to frame {ref_frame}")
    """
    if not config.enabled:
        LOGGER.debug("Alignment disabled by configuration")
        return None

    _require_mdanalysis("trajectory alignment")

    from MDAnalysis.analysis import align

    from polyzymd.analysis.core.centroid import find_centroid_frame

    n_frames = len(universe.trajectory)
    if stop_frame is None:
        stop_frame = n_frames

    ref_frame_idx: int | None = None

    if config.reference_mode == "centroid":
        # Find most populated state via K-Means clustering
        LOGGER.info(
            f"Aligning trajectory to centroid (selection: '{config.selection}', "
            f"frames {start_frame}-{stop_frame})"
        )

        ref_frame_idx = find_centroid_frame(
            universe,
            selection=config.centroid_selection,
            start_frame=start_frame,
            stop_frame=stop_frame,
            verbose=False,  # We handle our own logging
        )
        LOGGER.info(f"Found centroid at frame {ref_frame_idx} (0-indexed)")

        # Align to centroid frame
        aligner = align.AlignTraj(
            universe,
            universe,
            select=config.selection,
            ref_frame=ref_frame_idx,
            in_memory=True,
        ).run()
        del aligner

    elif config.reference_mode == "average":
        # Compute average structure and align to it
        LOGGER.info(f"Aligning trajectory to average structure (selection: '{config.selection}')")

        average = align.AverageStructure(
            universe,
            universe,
            select=config.selection,
            ref_frame=start_frame,  # Initial reference for iterative averaging
        ).run()
        ref_universe = average.results.universe

        aligner = align.AlignTraj(
            universe,
            ref_universe,
            select=config.selection,
            in_memory=True,
        ).run()
        del aligner
        del average

        ref_frame_idx = None  # Average is synthetic, not a real frame
        LOGGER.info("Aligned to average structure (synthetic reference)")

    elif config.reference_mode == "frame":
        # Use user-specified frame (convert from 1-indexed to 0-indexed)
        ref_frame_idx = config.reference_frame - 1  # type: ignore[operator]

        # Validate frame number
        if ref_frame_idx < 0 or ref_frame_idx >= n_frames:
            raise ValueError(
                f"reference_frame={config.reference_frame} (1-indexed) is out of range. "
                f"Valid range: 1 to {n_frames}"
            )

        LOGGER.info(
            f"Aligning trajectory to frame {config.reference_frame} "
            f"(selection: '{config.selection}')"
        )

        aligner = align.AlignTraj(
            universe,
            universe,
            select=config.selection,
            ref_frame=ref_frame_idx,
            in_memory=True,
        ).run()
        del aligner

    elif config.reference_mode == "external":
        # Align to an external PDB structure (e.g., catalytically competent crystal)
        import MDAnalysis as mda

        ref_path = Path(config.reference_file)  # type: ignore[arg-type]
        LOGGER.info(
            f"Aligning trajectory to external PDB: {ref_path} (selection: '{config.selection}')"
        )

        # Load external reference as a separate Universe
        ref_universe = mda.Universe(str(ref_path))

        # Validate atom counts match for the alignment selection
        traj_atoms = universe.select_atoms(config.selection)
        ref_atoms = ref_universe.select_atoms(config.selection)

        if len(ref_atoms) == 0:
            raise ValueError(
                f"External PDB '{ref_path.name}' has no atoms matching "
                f"selection '{config.selection}'. Check that the PDB contains "
                f"protein atoms with matching names/residues."
            )

        if len(traj_atoms) != len(ref_atoms):
            raise ValueError(
                f"Atom count mismatch between trajectory ({len(traj_atoms)}) "
                f"and external PDB ({len(ref_atoms)}) for selection "
                f"'{config.selection}'. The external PDB must contain "
                f"the same protein atoms as the simulation topology."
            )

        # Warn if residue names don't match
        traj_resnames = [r.resname for r in traj_atoms.residues]
        ref_resnames = [r.resname for r in ref_atoms.residues]
        if traj_resnames != ref_resnames:
            mismatches = [
                (i, t, r) for i, (t, r) in enumerate(zip(traj_resnames, ref_resnames)) if t != r
            ]
            LOGGER.warning(
                f"Residue name mismatches between trajectory and external PDB "
                f"({len(mismatches)} residues differ). First 5: "
                f"{mismatches[:5]}. Proceeding with alignment anyway."
            )

        LOGGER.info(
            f"External PDB matched: {len(ref_atoms)} atoms, {len(ref_atoms.residues)} residues"
        )

        # Align trajectory to external reference
        aligner = align.AlignTraj(
            universe,
            ref_universe,
            select=config.selection,
            in_memory=True,
        ).run()
        del aligner

        ref_frame_idx = None  # External reference is not a trajectory frame
        LOGGER.info("Aligned to external PDB reference structure")

    else:
        raise ValueError(f"Unknown reference_mode: {config.reference_mode}")

    LOGGER.info("Trajectory alignment complete")
    return ref_frame_idx


def get_alignment_description(config: AlignmentConfig) -> str:
    """Generate a human-readable description of the alignment configuration.

    Parameters
    ----------
    config : AlignmentConfig
        Alignment configuration.

    Returns
    -------
    str
        Description suitable for documentation or result files.
    """
    if not config.enabled:
        return "Alignment: disabled"

    if config.reference_mode == "centroid":
        return f"Alignment: centroid mode (K-Means clustering), selection='{config.selection}'"
    elif config.reference_mode == "average":
        return f"Alignment: average structure, selection='{config.selection}'"
    elif config.reference_mode == "frame":
        return (
            f"Alignment: frame {config.reference_frame} (1-indexed), selection='{config.selection}'"
        )
    elif config.reference_mode == "external":
        return f"Alignment: external PDB '{config.reference_file}', selection='{config.selection}'"
    else:
        return f"Alignment: {config.reference_mode}"
