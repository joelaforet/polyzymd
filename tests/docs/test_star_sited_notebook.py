"""The star-sited monomer notebook runs top to bottom.

``docs/source/tutorials/star_sited_monomers.ipynb`` is the sandbox the docs
and the polymer errors point users to. Its code cells run here in order, as a
script (the environment has no Jupyter kernel to drive), so a renamed import
or a changed builder API breaks this test instead of the user's first session.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

NOTEBOOK = (
    Path(__file__).resolve().parents[2]
    / "docs"
    / "source"
    / "tutorials"
    / "star_sited_monomers.ipynb"
)


def _code_cells() -> list[str]:
    """Return the source of each code cell of the notebook, in order."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]


def test_notebook_has_no_stored_outputs() -> None:
    """The notebook ships without outputs, so the docs show no stale results."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            assert cell["outputs"] == []
            assert cell["execution_count"] is None


@pytest.mark.skipif(shutil.which("packmol") is None, reason="needs the build environment")
def test_notebook_runs_top_to_bottom(capsys: pytest.CaptureFixture[str]) -> None:
    """Every cell runs, and the result is the star-sited SBMA, its end form and a built 3-mer."""
    pytest.importorskip("mbuild")
    namespace: dict[str, object] = {"__name__": "__notebook__"}
    for source in _code_cells():
        exec(compile(source, str(NOTEBOOK), "exec"), namespace)  # noqa: S102
    printed = capsys.readouterr().out

    sbma = "C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
    assert namespace["repeat_unit"] == f"*CC(*)(C){sbma}"
    assert namespace["end_cap"] == f"*C=C(C){sbma}"
    assert f'smiles: "*CC(*)(C){sbma}"' in printed
    assert "length: 5" in printed
    # An SBMA end, one SBMA repeat unit and an SBMA end: 38 + 39 + 38 atoms.
    assert "115 atoms, total charge +0.000 e" in printed


def _run_first_cells(smiles: str) -> dict[str, object]:
    """Run cells 1-3 of the notebook with ``smiles`` as the pasted monomer."""
    namespace: dict[str, object] = {"__name__": "__notebook__"}
    first, second, third = _code_cells()[:3]
    first = first.replace('"CC(=C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"', repr(smiles))
    for source in (first, second, third):
        exec(compile(source, str(NOTEBOOK), "exec"), namespace)  # noqa: S102
    return namespace


@pytest.mark.parametrize(
    ("smiles", "repeat_unit"),
    [
        ("C=Cc1c[nH]cn1", "*CC(*)c1c[nH]cn1"),
        ("CC(=C)C(=O)OC", "*CC(*)(C)C(=O)OC"),
        ("CC(C)NC(=O)C=C", "*CC(*)C(=O)NC(C)C"),
    ],
    ids=["aromatic_nH", "pubchem_order", "acrylamide"],
)
def test_the_default_linking_atoms_give_the_repeat_unit(smiles: str, repeat_unit: str) -> None:
    """The listed C=C, CH2 first, gives the repeat unit whatever the atom order."""
    assert _run_first_cells(smiles)["repeat_unit"] == repeat_unit
