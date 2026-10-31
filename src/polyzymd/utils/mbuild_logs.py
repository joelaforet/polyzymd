"""Mute mBuild's log while PolyzyMD replays modifications it does not simulate.

An unrelaxed replay (the audit's rebuild, ``validate --plans``) leaves the
close contacts mBuild warns about; those warnings describe the replay, not
the built structure, so they are kept out of the build log.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager


@contextmanager
def quiet_mbuild() -> Iterator[None]:
    """Raise every ``mbuild`` logger to ERROR for the duration, then restore them."""
    import mbuild.biopolymers  # noqa: F401  (registers mBuild's loggers)

    loggers = [
        logging.getLogger(name)
        for name in list(logging.root.manager.loggerDict)
        if name == "mbuild" or name.startswith("mbuild.")
    ]
    levels = [log.level for log in loggers]
    for log in loggers:
        log.setLevel(logging.ERROR)
    try:
        yield
    finally:
        for log, level in zip(loggers, levels):
            log.setLevel(level)
