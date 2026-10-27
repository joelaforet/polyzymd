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


def test_a_malformed_smiles_with_a_dot_is_a_value_error() -> None:
    """A SMILES RDKit cannot read names the co-solvent and the SMILES."""
    with pytest.raises(ValueError, match=r"surf.*CC\(=O\.\[Na\+\]"):
        solvent_molecules.split_counter_ions("CC(=O.[Na+]", name="surf")


def test_a_custom_smiles_under_a_library_name_gives_the_custom_molecule(empty_caches) -> None:
    """The bundled file is used only for the library molecule's own SMILES."""
    custom = solvent_molecules.get_solvent_molecule("dmso", smiles="CCO")
    assert custom.n_atoms == 9
    assert len(list(empty_caches.glob("dmso.*.sdf"))) == 1
    library = solvent_molecules.get_solvent_molecule("dmso", smiles="CS(C)=O")
    assert library is solvent_molecules.get_solvent_molecule("dmso")
    assert library.n_atoms == 10
    assert len(list(empty_caches.glob("*.sdf"))) == 1


def test_clear_cache_forgets_a_name_with_a_dot(empty_caches) -> None:
    """clear_cache(name) drops the in-memory molecule even when the name holds a dot."""
    solvent_molecules.get_solvent_molecule("my.solv", smiles="CCO")
    solvent_molecules.clear_cache("my.solv")
    assert solvent_molecules._loaded_molecules == {}
    assert list(empty_caches.glob("*.sdf")) == []
