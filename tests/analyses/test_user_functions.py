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


def test_same_named_files_of_two_studies_are_two_modules(tmp_path: Path) -> None:
    """Each study's f.py is its own module, so its dataclasses resolve their own names."""
    from polyzymd.analyses.user_functions import load_function

    code = (
        "from __future__ import annotations\n\nimport dataclasses\n\n\n"
        "@dataclasses.dataclass\nclass Value:\n    v: float\n\n\n"
        "def f():\n    return Value({v}).v\n"
    )
    functions = {}
    for study, value in (("lipa", 1.0), ("lipb", 2.0)):
        root = tmp_path / study
        (root / "analyses").mkdir(parents=True)
        (root / "study.yaml").write_text("equilibration: 0ns\n")
        (root / "analyses" / "f.py").write_text(code.format(v=value))
        functions[study] = load_function(root / "analyses" / "f.py", "f")
    assert functions["lipa"].__module__ == "polyzymd_study.lipa.f"
    assert functions["lipb"].__module__ == "polyzymd_study.lipb.f"
    assert functions["lipa"]() == 1.0 and functions["lipb"]() == 2.0
