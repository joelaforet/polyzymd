"""Load one replicate of a PolyzyMD run the way the analysis framework does.

:func:`load_replicate` returns the MDAnalysis universe and the production
window that every PolyzyMD analysis receives: restart segments joined into one
trajectory after their times are checked to line up, unfinished segments left
out, and the equilibration window resolved against the recorded frame times.
Use it to write a script, to try a measurement before turning it into an
analysis, or to call a plugin's ``compute()`` by hand::

    from polyzymd.analyses import iter_frames, load_replicate

    universe, frames = load_replicate("conditions/sbma/config.yaml", 1, equilibration="10ns")
    protein = universe.select_atoms("protein")
    rg = [protein.radius_of_gyration() for _ in iter_frames(universe, frames)]

``frames.run_kwargs()`` gives the same window as keyword arguments for an
MDAnalysis ``AnalysisBase.run()`` call.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple

if TYPE_CHECKING:
    from polyzymd.analyses.mda.frame_selection import FrameSelection


class Replicate(NamedTuple):
    """A loaded replicate and its production window."""

    universe: Any
    frames: FrameSelection


def load_replicate(
    config: Path | str | Any,
    replicate: int,
    *,
    equilibration: str,
    require_complete: bool = True,
    pbc_policy: str = "as_is",
) -> Replicate:
    """Load one replicate and resolve its production window.

    Parameters
    ----------
    config : Path, str or SimulationConfig
        The condition's ``config.yaml``, or the loaded config.
    replicate : int
        One-indexed replicate number.
    equilibration : str
        Time discarded from the start of the trajectory, for example
        ``"10ns"``. There is no default, so a script states its window.
    require_complete : bool, optional
        Leave out production segments the engine records as still running or
        failed, by default True.
    pbc_policy : str, optional
        ``"as_is"`` (default) reads coordinates as stored; ``"make_whole"``
        unwraps the protein and polymers and needs a topology with bonds.

    Returns
    -------
    Replicate
        ``(universe, frames)``, where ``frames`` is the
        :class:`~polyzymd.analyses.mda.frame_selection.FrameSelection` a
        plugin's ``compute()`` receives.

    Raises
    ------
    FileNotFoundError
        If the config, topology or trajectories are missing.
    TrajectoryLineageError
        If the restart segments overlap, run backwards or leave gaps.
    """
    from polyzymd.analyses.mda.frame_selection import FrameSelection
    from polyzymd.analyses.mda.universe import UniverseProvider
    from polyzymd.analyses.shared.loader import TrajectoryLoader
    from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

    if isinstance(config, (str, Path)):
        from polyzymd.config.loader import load_config

        config = load_config(config)
    loader = TrajectoryLoader(config)
    provider = UniverseProvider.from_config(
        config, loader=loader, require_complete=require_complete, pbc_policy=pbc_policy
    )
    universe = provider.load_universe(replicate)
    window = resolve_replicate_trajectory_window(
        loader=loader,
        replicate=replicate,
        equilibration=equilibration,
        n_frames_total=len(universe.trajectory),
    )
    return Replicate(universe, FrameSelection.from_trajectory_window(window))
