"""The replicate-loading notebook runs top to bottom in a new project.

``load_replicates.ipynb`` (shipped in ``polyzymd/templates/notebooks`` and
copied into each project by ``polyzymd project init``) is the sandbox of the
docs page "Analyze the trajectories with your own tools". Its code cells run
here in order, as a script, from the ``notebooks/`` folder of a project whose
study has one condition with a small OpenMM replicate, so a changed loader or
study API breaks this test instead of the user's first session.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.project_scaffold import PROJECT_NOTEBOOKS

NOTEBOOKS = Path(__file__).resolve().parents[2] / "src" / "polyzymd" / "templates" / "notebooks"
NOTEBOOK = NOTEBOOKS / "load_replicates.ipynb"

#: Two alanine backbones, N CA C O each, in Å.
BACKBONE = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.0, 1.4, 0.0],
        [1.4, 2.4, 0.0],
        [3.3, 1.5, 0.0],
        [4.0, 2.7, 0.0],
        [5.5, 2.5, 0.0],
        [6.1, 1.4, 0.0],
    ]
)


def _code_cells() -> list[str]:
    """Return the source of each code cell of the notebook, in order."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]


def _run(namespace: dict[str, object]) -> None:
    for source in _code_cells():
        exec(compile(source, str(NOTEBOOK), "exec"), namespace)  # noqa: S102


def _project_with_a_replicate(root: Path) -> Path:
    """Write a project whose study ``trpcage`` has the condition ``Water`` with replicate 1."""
    from polyzymd.analyses.project_scaffold import create_project
    from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

    created = create_project(root, ["trpcage"], git=False)
    folder = created.root / "trpcage" / "conditions" / "water"
    config = write_simulation_config(folder, scratch=created.root / "runs" / "water")
    (folder / "test.pdb").write_text("REMARK input\nEND\n")
    write_openmm_frames(
        config,
        1,
        [BACKBONE * scale for scale in (1.0, 1.1, 1.2, 1.3, 1.4)],
        [0, 0, 0, 0, 1, 1, 1, 1],
        names=["N", "CA", "C", "O"] * 2,
        resnames=["ALA", "ALA"],
        elements=["N", "C", "C", "O"] * 2,
        chain_ids=["A"] * 8,
    )
    study_yaml = created.root / "trpcage" / "study.yaml"
    text = study_yaml.read_text().replace(
        "conditions: {}", 'conditions:\n  "Water": "conditions/water/config.yaml"'
    )
    study_yaml.write_text(text)
    return created.root


def test_the_notebook_loads_the_first_replicate_of_the_project(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Run from notebooks/, it loads Water replicate 1 both ways and measures it."""
    pytest.importorskip("MDAnalysis")
    pytest.importorskip("matplotlib").use("Agg")
    root = _project_with_a_replicate(tmp_path / "paper")
    monkeypatch.chdir(root / "notebooks")

    namespace: dict[str, object] = {"__name__": "__notebook__"}
    _run(namespace)
    printed = capsys.readouterr().out

    assert "Chosen: Replicate(condition='Water', index=1)" in printed
    assert "5 production frames; 5 after the 0ns window" in printed
    assert "protein                'protein'                              8 atoms" in printed
    assert "Same atoms as PolyzyMD's universe: True" in printed
    assert "Same first frame: True" in printed
    by_hand = namespace["by_hand"]
    expected = [by_hand.atoms.radius_of_gyration() for _ in by_hand.trajectory]
    np.testing.assert_allclose(namespace["values"], expected, rtol=1e-5)


def test_a_new_project_is_told_to_add_a_condition(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """In a project fresh from init, the first cell says what to do next."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project_scaffold import create_project

    created = create_project(tmp_path / "paper", ["trpcage"], git=False)
    monkeypatch.chdir(created.root / "notebooks")
    with pytest.raises(ProtocolError, match="polyzymd study add-condition"):
        exec(compile(_code_cells()[0], str(NOTEBOOK), "exec"), {})  # noqa: S102


@pytest.mark.parametrize("name", PROJECT_NOTEBOOKS)
def test_each_shipped_notebook_has_no_outputs_and_reaches_new_projects(
    name: str, tmp_path: Path
) -> None:
    """Every notebook ships without outputs, and project init copies it byte for byte."""
    from polyzymd.analyses.project_scaffold import create_project

    shipped = NOTEBOOKS / name
    for cell in json.loads(shipped.read_text(encoding="utf-8"))["cells"]:
        if cell["cell_type"] == "code":
            assert cell["outputs"] == []
            assert cell["execution_count"] is None
    created = create_project(tmp_path / "paper", ["study_a"], git=False)
    assert (created.root / "notebooks" / name).read_bytes() == shipped.read_bytes()
