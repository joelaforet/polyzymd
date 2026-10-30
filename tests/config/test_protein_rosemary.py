"""A protein may opt in to OpenFF Rosemary alpha for the whole component.

Which charges a protein takes is read from its force field file's content
(protein library charges, or a NAGLCharges section), not from the file name.
"""

from __future__ import annotations

import copy
import shutil
import sys
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from polyzymd.config.schema import (
    ProteinSolute,
    SimulationConfig,
    is_rosemary,
    protein_charge_handlers,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.rosemary import ROSEMARY_ALPHA0, newest_rosemary  # noqa: E402

ROSEMARY = newest_rosemary()
NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}
QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"
ACYL = [
    {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
    {
        "attach": {
            "type": "generic",
            "site": "A:K8",
            "atom_name": "NZ",
            "fragment": {"smiles": "*C(=O)CCC", "resname": "BUT"},
        }
    },
]


def _protein(**keys) -> ProteinSolute:
    return ProteinSolute(type="protein", pdb_path="trpcage.pdb", **keys)


def _rosemary_copy(tmp_path: Path, name: str) -> str:
    """Rosemary alpha's file under another name or path."""
    from openforcefields import get_forcefield_dirs_paths

    source = next(
        Path(directory) / ROSEMARY
        for directory in get_forcefield_dirs_paths()
        if (Path(directory) / ROSEMARY).is_file()
    )
    target = tmp_path / name
    shutil.copy(source, target)
    return str(target)


def test_the_force_field_content_says_how_it_charges_a_protein() -> None:
    assert protein_charge_handlers(("ff14sb_off_impropers_0.0.4.offxml",)) == (None, True)
    for rosemary in {ROSEMARY, ROSEMARY_ALPHA0}:
        model, library = protein_charge_handlers((rosemary,))
        assert model is not None and not library
        assert is_rosemary([rosemary])
    assert not is_rosemary(["ff14sb_off_impropers_0.0.4.offxml"])


def test_the_first_rosemary_candidate_stays_accepted() -> None:
    assert _protein(force_field=[ROSEMARY_ALPHA0], charges=NAGL).on_rosemary


def test_a_rosemary_protein_takes_nagl_charges_and_its_modifications_need_no_parameters() -> None:
    protein = _protein(force_field=[ROSEMARY], charges=NAGL, modifications=ACYL)
    assert protein.on_rosemary and protein.is_modified


def test_a_rosemary_protein_refuses_library_charges_and_names_its_model() -> None:
    with pytest.raises(ValidationError, match="has no protein library charges; it charges with"):
        _protein(force_field=[ROSEMARY], charges={"method": "library"})


def test_a_future_named_rosemary_file_with_library_charges_is_refused(tmp_path) -> None:
    future = _rosemary_copy(tmp_path, "openff_no_water-3.0.0.offxml")
    with pytest.raises(ValidationError, match="Use charges: {method: nagl"):
        _protein(force_field=[future], charges={"method": "library"})


def test_a_custom_path_rosemary_file_with_nagl_charges_is_accepted(tmp_path) -> None:
    custom = _rosemary_copy(tmp_path, "my_rosemary.offxml")
    assert _protein(force_field=[custom], charges=NAGL).on_rosemary


def test_an_ff14sb_protein_keeps_library_charges_and_refuses_nagl() -> None:
    assert not _protein(
        force_field=["ff14sb_off_impropers_0.0.4.offxml"], charges={"method": "library"}
    ).on_rosemary
    with pytest.raises(ValidationError, match="has protein library charges"):
        _protein(force_field=["ff14sb_off_impropers_0.0.4.offxml"], charges=NAGL)


def test_a_rosemary_protein_refuses_modification_parameters() -> None:
    with pytest.raises(ValidationError, match="modification_parameters does not apply"):
        _protein(
            force_field=[ROSEMARY],
            charges=NAGL,
            modifications=ACYL,
            modification_parameters={
                "generic": {"force_field": ["openff-2.3.0.offxml"], "charges": NAGL}
            },
        )


def test_rosemary_is_the_protein_force_field_on_its_own() -> None:
    with pytest.raises(ValidationError, match="only force field file"):
        _protein(force_field=[ROSEMARY, "openff-2.3.0.offxml"], charges=NAGL)


def test_the_protein_and_a_co_solvent_share_one_nagl_model() -> None:
    data = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    data["solutes"][0].update(
        pdb_path=str(QUICKSTART / "trpcage.pdb"), force_field=[ROSEMARY], charges=NAGL
    )
    data["solvent"]["co_solvents"] = [
        {"name": "x", "smiles": "CCO", "count": 1, "charges": {**NAGL, "model": "b.pt"}}
    ]
    with pytest.raises(ValidationError, match="whole protein on its NAGL force field") as error:
        SimulationConfig.model_validate(copy.deepcopy(data))
    assert "one model per build" in str(error.value)
