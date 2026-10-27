"""Replicates run with the user's own protocol, outside PolyzyMD.

A user may build a system with PolyzyMD and simulate it any way they like (an
OpenMM script with PLUMED, their own GROMACS ``.mdp`` files, enhanced
sampling), then put the trajectory back into the replicate folder for
analysis. They declare it by copying the inputs they ran into ``protocol/``:
the replicate folder's for OpenMM, ``gromacs/protocol/`` for GROMACS.

PolyzyMD then vouches for the system it built (its audited build record),
not for the run. It never runs, resubmits or recovers such a replicate, so
the user's trajectories are never overwritten or deleted; ``polyzymd status``
shows it as run with its own protocol; analysis records say so; and a study
deposit holds the protocol files, never the run parameters PolyzyMD wrote
but did not run.
"""

from __future__ import annotations

import os
from pathlib import Path

from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

#: The folder that holds the inputs of a run made outside PolyzyMD.
PROTOCOL_DIR = "protocol"


class OwnProtocolError(ArtifactIntegrityError):
    """Raised when PolyzyMD is asked to run a replicate run with the user's own protocol."""


def protocol_dir(working_dir: str | Path) -> Path | None:
    """Return the replicate's ``protocol/`` folder, if it holds any file.

    ``working_dir`` is the replicate folder or its ``gromacs/`` folder; both
    ``protocol/`` and ``gromacs/protocol/`` are found from either.
    """
    folder = Path(working_dir)
    replicate = folder.parent if folder.name == "gromacs" else folder
    for candidate in (replicate / PROTOCOL_DIR, replicate / "gromacs" / PROTOCOL_DIR):
        if candidate.is_dir() and any(path.is_file() for path in candidate.rglob("*")):
            return candidate
    return None


def protocol_files(working_dir: str | Path) -> list[Path]:
    """Return every file under the replicate's ``protocol/`` folder, sorted."""
    folder = protocol_dir(working_dir)
    if folder is None:
        return []
    return sorted(path for path in folder.rglob("*") if path.is_file())


def require_polyzymd_run(working_dir: str | Path, action: str) -> None:
    """Refuse to ``action`` a replicate that was run with the user's own protocol."""
    folder = protocol_dir(working_dir)
    if folder is None:
        return
    folder = Path(os.path.normpath(folder))
    raise OwnProtocolError(
        f"{folder.parent} was run with your own protocol ({folder}): PolyzyMD does not "
        f"{action} it, so its trajectories are never overwritten or deleted.\n"
        "fix: analyze it with polyzymd analyze. To run this replicate with PolyzyMD "
        "instead, move protocol/ and the trajectories out of the replicate folder first."
    )
