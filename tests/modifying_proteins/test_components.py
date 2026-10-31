"""The component record of a modified protein."""

from pathlib import Path

from polyzymd import modifying_proteins as mp
from polyzymd.modifying_proteins.components import protein_component

TRPCAGE = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"
SAGE = mp.ParameterSet(
    force_field=("openff-2.3.0.offxml",),
    charges={"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
)


def test_the_protein_keeps_the_segment_id_its_config_gives() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    mp.attach(
        protein,
        {"smiles": "*C(=O)C", "resname": "ACE"},
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )
    record = mp.parameterize(protein, {"generic": SAGE}).record
    default = protein_component(protein, record, {"generic": SAGE})
    custom = protein_component(protein, record, {"generic": SAGE}, segid="ENZ")
    assert (default.segid, custom.segid) == ("PROT", "ENZ")
    assert set(custom.segments) == {"M001"}


def test_the_rebuild_recipe_names_inputs_relative_to_the_config_folder(tmp_path) -> None:
    import hashlib
    import shutil

    folder = tmp_path / "study"
    (folder / "structures").mkdir(parents=True)
    cif = shutil.copy(
        Path(__file__).resolve().parents[1] / "data" / "modifying_proteins" / "4II.cif",
        folder / "structures",
    )
    shutil.copy(TRPCAGE, folder)
    protein = mp.load_protein(folder / "trpcage.pdb")
    mp.mutate(protein, "A:Y3", {"ccd_file": str(cif)}, relax=False)
    record = mp.parameterize(protein, {"noncanonical": SAGE}).record
    component = protein_component(protein, record, {"noncanonical": SAGE}, base_dir=folder)
    recipe = component.domains[0].reference["modified_protein"]
    (mutation,) = recipe["spec"]["mutations"]
    assert mutation["to"]["ccd_file"] == "structures/4II.cif"
    digest = hashlib.sha256(Path(cif).read_bytes()).hexdigest()
    assert recipe["files"]["structures/4II.cif"] == digest
    assert "trpcage.pdb" in recipe["files"]


def test_identical_chains_share_one_charge_template(tmp_path) -> None:
    """A homodimer with the same modification on both chains parameterizes (Interchange
    refuses two isomorphic charge templates)."""
    atoms = [
        line for line in TRPCAGE.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    chain_b = [
        f"{line[:21]}B{line[22:30]}{float(line[30:38]) + 30:8.3f}{line[38:]}" for line in atoms
    ]
    pdb = tmp_path / "dimer.pdb"
    pdb.write_text("\n".join(atoms + ["TER"] + chain_b + ["TER", "END"]) + "\n")
    protein = mp.load_protein(pdb)
    for chain in "AB":
        mp.deprotonate(protein, f"{chain}:K8", "NZ")
        mp.attach(
            protein,
            {"smiles": "*C(=O)C", "resname": "ACE"},
            type="generic",
            site=f"{chain}:K8",
            atom_name="NZ",
            relax=False,
        )
    result = mp.parameterize(protein, {"generic": SAGE})
    assert len(result.record.molecules) == 2
    assert result.interchange.topology.n_molecules == 2
