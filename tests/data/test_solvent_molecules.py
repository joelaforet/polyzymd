"""Tests for co-solvent molecules made from SMILES and their disk cache."""

import pytest

from polyzymd.data import solvent_molecules


@pytest.fixture()
def empty_caches(tmp_path, monkeypatch):
    """Point the co-solvent cache at an empty folder and empty the in-memory cache."""
    monkeypatch.setattr(solvent_molecules, "_USER_CACHE_DIR", tmp_path)
    monkeypatch.setattr(solvent_molecules, "_loaded_molecules", {})
    return tmp_path


def test_a_changed_smiles_with_the_same_name_gives_the_new_molecule(empty_caches, monkeypatch):
    """The cache follows the SMILES: a new SMILES under an old name is not served the old molecule."""
    first = solvent_molecules.get_solvent_molecule("surf", smiles="CCO")
    assert first.n_atoms == 9
    # A new session: only the disk cache is left.
    monkeypatch.setattr(solvent_molecules, "_loaded_molecules", {})
    second = solvent_molecules.get_solvent_molecule("surf", smiles="CCCO")
    assert second.n_atoms == 12
    assert solvent_molecules.get_solvent_molecule("surf", smiles="CCCO") is second
    assert len(list(empty_caches.glob("*.sdf"))) == 2


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCCCCCCCCCCCOS(=O)(=O)[O-]", ("CCCCCCCCCCCCOS(=O)(=O)[O-]", 0, 0)),
        ("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", ("CCCCCCCCCCCCOS(=O)(=O)[O-]", 1, 0)),
        ("[Na+].[Na+].O=C([O-])CCC(=O)[O-]", ("O=C([O-])CCC(=O)[O-]", 2, 0)),
        ("C[N+](C)(C)C.[Cl-]", ("C[N+](C)(C)C", 0, 1)),
    ],
)
def test_counter_ions_are_split_off_a_smiles(smiles, expected) -> None:
    """Na+ and Cl- written in a SMILES are counted and removed from the molecule."""
    assert solvent_molecules.split_counter_ions(smiles) == expected


@pytest.mark.parametrize("smiles", ["CCCCCCCCCCCCOS(=O)(=O)[O-].[K+]", "CCO.CO"])
def test_other_parts_of_a_smiles_are_refused(smiles) -> None:
    """A SMILES with a part that is not the molecule, Na+ or Cl- is refused with a hint."""
    with pytest.raises(ValueError, match="own co_solvents entry"):
        solvent_molecules.split_counter_ions(smiles)
