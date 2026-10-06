"""Study, Condition and Replicate: how a replicate reports a universe it cannot load."""

from __future__ import annotations

from pathlib import Path

import pytest


def test_a_topology_that_cannot_load_is_named_as_such(tmp_path: Path) -> None:
    """A broken topology is not reported as an equilibration problem."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study import Condition, Replicate

    class Broken:
        def load_universe(self, index):
            raise ValueError("Length of charges does not match number of atoms")

    replicate = Replicate.__new__(Replicate)
    replicate._universe = None
    replicate.index = 1
    replicate.condition = type("C", (), {"_provider": Broken(), "label": "A"})()
    with pytest.raises(ProtocolError, match="Cannot load the topology") as info:
        Replicate.universe(replicate)
    assert "analysis-topology" in info.value.hint
    assert Condition  # the class used above is the study's own
