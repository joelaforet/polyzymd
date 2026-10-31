"""One force field for the whole protein (OpenFF Rosemary), NAGL per molecule."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.rosemary import newest_rosemary  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}

pytestmark = pytest.mark.filterwarnings("ignore")


def _two_chains(path: Path) -> Path:
    """Trp-cage as chain A and, 3 nm away, chain B."""
    atoms = [
        line for line in TRPCAGE.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    shifted = [
        f"{line[:21]}B{line[22:30]}{float(line[30:38]) + 30:8.3f}{line[38:]}" for line in atoms
    ]
    path.write_text("\n".join(atoms + ["TER"] + shifted + ["TER", "END"]) + "\n")
    return path


def test_each_chain_is_charged_to_its_own_formal_charge(tmp_path) -> None:
    """NAGL runs per molecule, so a chain's charges sum to that chain's formal charge."""
    from openff.units import unit

    from polyzymd import modifying_proteins as mp

    protein = mp.load_protein(_two_chains(tmp_path / "two.pdb"))
    mp.deprotonate(protein, "A:K8", "NZ")  # chain A: +1 -> 0; chain B stays +1
    result = mp.parameterize_whole(protein, (newest_rosemary(),), NAGL)
    first, second = result.record.molecules
    charges = result.interchange["Electrostatics"].get_charge_array().m_as(unit.elementary_charge)
    assert np.isclose(charges[first[0] : first[1]].sum(), 0.0, atol=1e-4)
    assert np.isclose(charges[second[0] : second[1]].sum(), 1.0, atol=1e-4)
    assert result.interchange.topology.n_molecules == 2
    assert result.record.protein_charges == NAGL


def _neutral_n_terminus(protein) -> None:
    from polyzymd import modifying_proteins as mp

    mp.deprotonate(protein, "A:N1", "N")


def _protonated_c_terminus(protein) -> None:
    from polyzymd import modifying_proteins as mp

    mp.protonate(protein, "A:S20", "OXT")


@pytest.mark.parametrize("change", [_neutral_n_terminus, _protonated_c_terminus])
def test_a_terminus_changed_by_protonation_alone_is_parameterized_whole(change) -> None:
    """Rosemary needs no residue templates, so a neutral terminus builds."""
    from openff.units import unit

    from polyzymd import modifying_proteins as mp

    protein = mp.load_protein(TRPCAGE)
    change(protein)
    result = mp.parameterize_whole(protein, (newest_rosemary(),), NAGL)
    charges = result.interchange["Electrostatics"].get_charge_array().m_as(unit.elementary_charge)
    # Trp-cage is +1 with charged termini; a neutral terminus removes one charge.
    expected = 0.0 if change is _neutral_n_terminus else 2.0
    assert np.isclose(charges.sum(), expected, atol=1e-4)


def test_a_labelled_terminus_is_parameterized_whole() -> None:
    from polyzymd import modifying_proteins as mp

    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:N1", "N")
    mp.attach(
        protein,
        {"smiles": "*C(C)=O", "resname": "ACE"},
        type="generic",
        site="A:N1",
        atom_name="N",
        relax=False,
    )
    result = mp.parameterize_whole(protein, (newest_rosemary(),), NAGL)
    assert result.record.protein_charges == NAGL
