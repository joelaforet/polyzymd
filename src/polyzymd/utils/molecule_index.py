"""Constant-time ``Topology.molecule_index`` while exporting large systems.

The OpenFF toolkit's ``Topology.molecule_index`` scans the topology's molecules
on every call. OpenFF Interchange calls it once per molecule when it writes an
OpenMM topology and the nonbonded part of an OpenMM system, so exporting a
solvated system costs time that grows with the square of the number of
molecules: about ten minutes for 106,000 molecules (a glycosylated BG505
protomer in water), against about thirty seconds when the lookup is constant
time.

:func:`indexed_molecules` replaces the method, for the duration of a ``with``
block, with one that looks the molecule up in an identity map built once per
topology. Every index the map gives is checked against the topology
(``topology.molecule(index) is molecule``, constant time) before it is
returned; on a mismatch the map is rebuilt once and checked again, and a
molecule still not found goes to the toolkit's own scan. The index is therefore
the one the toolkit gives for the topology as it is at the time of the call,
including after molecules are added, removed or replaced inside the block, and
a molecule that is not in the topology raises the toolkit's error. Remove this
module when the toolkit's lookup is constant time (tracked in issue #281).

The patch is process-wide. Overlapping or nested blocks share one patch and a
depth counter, so the toolkit's method is restored when the last block exits,
whatever the order the blocks exit in. A lock guards the counter, but the patch
is not otherwise thread-safe: a thread outside any block still sees the patched
method while another thread is inside one. PolyzyMD builds and exports run in
one thread.
"""

from __future__ import annotations

import contextlib
import threading
from typing import Any, Iterator

_lock = threading.Lock()
_depth = 0
_original: Any = None
# id(topology) -> (topology, {id(molecule): index}); the topology is kept so a
# reused id cannot select another topology's map.
_maps: dict[int, tuple[Any, dict[int, int]]] = {}


def _build(topology: Any) -> dict[int, int]:
    index = {id(m): i for i, m in enumerate(topology.molecules)}
    _maps[id(topology)] = (topology, index)
    return index


def _holds(topology: Any, index: int | None, molecule: Any) -> bool:
    return (
        index is not None and index < topology.n_molecules and topology.molecule(index) is molecule
    )


def _molecule_index(self: Any, molecule: Any) -> int:
    entry = _maps.get(id(self))
    index_map = entry[1] if entry is not None and entry[0] is self else _build(self)
    index = index_map.get(id(molecule))
    if _holds(self, index, molecule):
        return index
    # The topology changed since the map was built, or the molecule is not in it.
    index = _build(self).get(id(molecule))
    if _holds(self, index, molecule):
        return index
    return _original(self, molecule)


@contextlib.contextmanager
def indexed_molecules() -> Iterator[None]:
    """Make ``Topology.molecule_index`` constant time inside the ``with`` block.

    Inside the block the method returns the index the toolkit's own method gives
    for the topology as it is at the time of the call, including after the
    topology is changed inside the block, and raises the toolkit's
    ``MoleculeNotInTopologyError`` for a molecule that is not in it. Blocks may
    overlap; the toolkit's method is restored when the last one exits. Not
    thread-safe beyond that (see the module docstring).
    """
    global _depth, _original
    from openff.toolkit import Topology

    with _lock:
        if _depth == 0:
            _original = Topology.molecule_index
            Topology.molecule_index = _molecule_index
        _depth += 1
    try:
        yield
    finally:
        with _lock:
            _depth -= 1
            if _depth == 0:
                Topology.molecule_index = _original
                _original = None
                _maps.clear()
