"""Test an analysis on in-memory universes, without trajectories or disk.

:func:`run_in_memory` runs a plugin's ``compute()`` once per replicate universe
and reduces and aggregates the result with the same functions the framework
uses, so a known-answer test checks both the measurement and the statistics
the analysis will report. :func:`synthetic_universe` builds a universe whose
answer is known::

    from polyzymd.analyses.testing import run_in_memory, synthetic_universe

    aggregates = run_in_memory(MyAnalysis, {"selection": "all"}, synthetic_universe())
    assert aggregates[0].mean == pytest.approx(1.0)

``polyzymd new-analysis`` writes tests in this form.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Sequence

if TYPE_CHECKING:
    from polyzymd.analyses.contract import ObservableAggregate
    from polyzymd.analyses.mda.frame_selection import FrameSelection

#: Four atoms on a unit cross. The radius of gyration of this shape is its scale.
CROSS = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))


def synthetic_universe(scale: float = 1.0, n_frames: int = 5) -> Any:
    """Four unit-mass atoms on a cross, identical on every frame.

    Parameters
    ----------
    scale : float, optional
        Distance of each atom from the origin, by default 1.0. The radius of
        gyration of the group is exactly this value.
    n_frames : int, optional
        Number of frames, by default 5.

    Returns
    -------
    MDAnalysis.Universe
        Universe backed by ``MemoryReader``.
    """
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0] * 4)
    positions = np.asarray(CROSS, dtype=np.float32) * scale
    universe.load_new(np.stack([positions] * n_frames), format=MemoryReader)
    return universe


def all_frames(universe: Any) -> FrameSelection:
    """A production window covering every frame of ``universe``."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    n_frames = len(universe.trajectory)
    return FrameSelection(start=0, stop=n_frames, step=1, n_frames_total=n_frames)


def run_in_memory(
    analysis: Any,
    settings: Any,
    universes: Any,
    *,
    n_replicates: int = 3,
) -> list[ObservableAggregate]:
    """Run an analysis on in-memory replicates and aggregate them.

    Parameters
    ----------
    analysis : plugin class, plugin instance or type[Analysis]
        The analysis under test, as a study module defines it or as
        :func:`~polyzymd.analyses.contract.contract_analysis` returned it.
    settings : BaseModel or dict
        Settings, validated against the plugin's ``Settings`` model.
    universes : Universe or sequence of Universe
        One universe per replicate, or one universe used for every replicate.
    n_replicates : int, optional
        Replicates to run when a single universe is given, by default 3.

    Returns
    -------
    list[ObservableAggregate]
        One aggregate per observable, as the framework reports it for one
        condition.

    Raises
    ------
    PluginContractError
        If ``compute()`` returns something other than observables, or the
        replicates disagree on what they report.
    """
    from polyzymd.analyses.base import _unpack
    from polyzymd.analyses.contract import aggregate_observables

    plugin = getattr(analysis, "plugin", None) or (
        analysis() if isinstance(analysis, type) else analysis
    )
    if not isinstance(settings, plugin.Settings):
        settings = plugin.Settings.model_validate(settings or {})
    if isinstance(universes, Sequence):
        replicates = list(universes)
    else:
        replicates = [universes] * n_replicates
    measured = [
        _unpack(plugin.compute(universe, all_frames(universe), settings))[0]
        for universe in replicates
    ]
    return aggregate_observables(measured)
