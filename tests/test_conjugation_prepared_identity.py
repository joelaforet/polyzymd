from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records,
    parse_pdb_conect_pairs,
)
from polyzymd.builders.conjugation.structure.prepared_identity import (
    canonicalize_prepared_protein_identity,
)


def test_prepared_identity_maps_multichain_atoms_and_connectivity(tmp_path) -> None:
    """Canonical output should retain boundaries, links, and source provenance."""
    from openmm.app import PDBFile

    source = tmp_path / "source.pdb"
    prepared = tmp_path / "prepared.pdb"
    output = tmp_path / "canonical.pdb"
    source.write_text(
        "".join(
            (
                _atom(1, "N", "A", 1, "N"),
                _atom(2, "C", "A", 1),
                _atom(3, "SG", "A", 1, "S"),
                "TER\n",
                _atom(4, "N", "B", 1, "N"),
                _atom(5, "C", "B", 1),
                _atom(6, "SG", "B", 1, "S"),
                _atom(7, "CB", "B", 1),
                "CONECT    3    6\nEND\n",
            )
        )
    )
    prepared.write_text(
        "".join(
            (
                _atom(10, "N", "A", 1, "N"),
                _atom(11, "H", "A", 1, "H"),
                _atom(12, "C", "A", 1),
                _atom(13, "SG", "A", 1, "S"),
                "TER\n",
                _atom(20, "N", "B", 1, "N"),
                _atom(21, "H", "B", 1, "H"),
                _atom(22, "C", "B", 1),
                _atom(23, "SG", "B", 1, "S"),
                "CONECT   10   12\nEND\n",
            )
        )
    )

    result = canonicalize_prepared_protein_identity(source, prepared, output)
    atoms = parse_pdb_atom_records(output)
    mapped = result.canonical_for_source(
        SimpleNamespace(
            chain_id="B",
            residue_name="CYS",
            residue_number=1,
            insertion_code="",
            atom_name="SG",
            atom_serial=6,
            atom_index=5,
        )
    )

    assert {(atom.chain_id, atom.residue_number) for atom in atoms} == {("A", 1), ("A", 2)}
    assert (mapped.chain_id, mapped.residue_number, mapped.serial, mapped.atom_index) == (
        "A",
        2,
        8,
        7,
    )
    assert sum(record.status == "prepared_added" for record in result.records) == 2
    assert result.mapping_path.exists()
    assert result.canonical_conect_pair_count == 2
    assert {tuple(sorted(pair)) for pair in parse_pdb_conect_pairs(output)} == {(1, 3), (4, 8)}
    assert output.read_text().count("TER\n") == 2
    topology = PDBFile(str(output)).topology
    residues = list(topology.residues())
    bonds = {frozenset(pair) for pair in topology.bonds()}
    assert frozenset((_named(residues[0], "C"), _named(residues[1], "N"))) not in bonds
    assert frozenset((_named(residues[0], "SG"), _named(residues[1], "SG"))) in bonds
    with pytest.raises(ValueError, match="did not survive"):
        result.canonical_for_source(
            SimpleNamespace(
                chain_id="B",
                residue_name="CYS",
                residue_number=1,
                insertion_code="",
                atom_name="CB",
                atom_serial=7,
                atom_index=6,
            )
        )

    repeated = tmp_path / "canonical_repeated.pdb"
    canonicalize_prepared_protein_identity(output, output, repeated)
    assert repeated.read_bytes() == output.read_bytes()


def test_prepared_identity_rejects_ambiguous_source_atom(tmp_path) -> None:
    """Repeated source atom identities must fail rather than map arbitrarily."""
    source = tmp_path / "ambiguous.pdb"
    source.write_text(_atom(1, "N", "A", 1) + _atom(2, "N", "A", 1) + "END\n")
    with pytest.raises(ValueError, match="Ambiguous source PDB atom identity"):
        canonicalize_prepared_protein_identity(source, source, tmp_path / "output.pdb")


def _atom(serial: int, name: str, chain: str, residue: int, element: str = "C") -> str:
    return (
        f"ATOM  {serial:5d} {name:>4} CYS {chain}{residue:4d}    "
        f"{serial:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2}\n"
    )


def _named(residue, name: str):
    return next(atom for atom in residue.atoms() if atom.name == name)
