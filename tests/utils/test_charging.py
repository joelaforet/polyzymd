"""Tests for the partial-charge methods."""

import pytest

from polyzymd.utils import charging


def test_am1bcc_without_ambertools_fails_with_the_fix(monkeypatch) -> None:
    """AM1-BCC without AmberTools stops with an error that names the fix, not an unexpected error."""
    from openff.toolkit import Molecule

    monkeypatch.setattr(charging, "am1bcc_available", lambda: False)
    with pytest.raises(ValueError, match="AmberTools.*nagl"):
        charging.get_charger("am1bcc").charge_molecule(Molecule.from_smiles("CCO"))


def test_nagl_charges_with_a_home_that_has_no_cache_folder_yet(monkeypatch, tmp_path) -> None:
    """A fresh home (macOS ~/Library/Caches missing) still charges with the installed model."""
    from openff.nagl_models import _dynamic_fetch
    from openff.toolkit import Molecule

    cache = tmp_path / "home" / "Library" / "Caches" / "OPENFF_NAGL_MODELS"
    monkeypatch.setattr(_dynamic_fetch, "CACHE_DIR", cache)
    molecule = charging.get_charger("nagl").charge_molecule(Molecule.from_smiles("CCO"))
    assert molecule.partial_charges is not None
    assert cache.is_dir()
