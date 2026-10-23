"""Tests for the partial-charge methods."""

import pytest

from polyzymd.utils import charging


def test_am1bcc_without_ambertools_fails_with_the_fix(monkeypatch) -> None:
    """AM1-BCC without AmberTools stops with an error that names the fix, not an unexpected error."""
    from openff.toolkit import Molecule

    monkeypatch.setattr(charging, "am1bcc_available", lambda: False)
    with pytest.raises(ValueError, match="AmberTools.*nagl"):
        charging.get_charger("am1bcc").charge_molecule(Molecule.from_smiles("CCO"))
