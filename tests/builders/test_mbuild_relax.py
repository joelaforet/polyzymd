"""The builder environments can relax an mBuild attachment."""

from pathlib import Path

TRPCAGE = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"


def test_attach_with_relax_runs_in_the_builder_environment() -> None:
    from mbuild.biopolymers import Protein, fragment_from_smiles

    protein = Protein(str(TRPCAGE))
    n_before = protein.n_particles
    protein.deprotonate(8, "NZ")
    fragment = fragment_from_smiles("*C(=O)CCCCCCC", "OC8")
    protein.attach(fragment, resnum=8, atom_name="NZ", relax=True)
    # Octanoyl adds C8H15O; deprotonation and the amide bond remove two H.
    assert protein.n_particles == n_before + 24 - 2
