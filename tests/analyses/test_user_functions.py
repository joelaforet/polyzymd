"""Loading a study's own function files."""

from __future__ import annotations

import pickle
from pathlib import Path

from polyzymd.analyses.user_functions import load_module


def test_a_dataclass_with_future_annotations_loads_and_pickles(tmp_path: Path) -> None:
    """dataclasses and pickle look the module up in sys.modules, so it is there."""
    file = tmp_path / "shapes.py"
    file.write_text(
        "from __future__ import annotations\n"
        "from dataclasses import dataclass\n\n\n"
        "@dataclass\n"
        "class Box:\n"
        "    side: float\n"
    )
    module = load_module(file)
    box = pickle.loads(pickle.dumps(module.Box(2.0)))
    assert box.side == 2.0
