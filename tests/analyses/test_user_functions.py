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


def test_each_study_imports_its_own_helper_module(tmp_path: Path) -> None:
    """Two studies' analyses/ folders with a helper.py each: each function sees its own."""
    from polyzymd.analyses.user_functions import load_function

    values = {}
    for study, value in (("lipa", 1), ("lipb", 2)):
        folder = tmp_path / study / "analyses"
        folder.mkdir(parents=True)
        (folder / "helper.py").write_text(f"V = {value}\n")
        (folder / "metric.py").write_text("import helper\n\n\ndef v():\n    return helper.V\n")
        values[study] = load_function(folder / "metric.py", "v")
    assert {study: function() for study, function in values.items()} == {"lipa": 1, "lipb": 2}
