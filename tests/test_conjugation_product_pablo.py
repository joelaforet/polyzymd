"""Tests for product-state Pablo residue-library construction."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from polyzymd.builders.conjugation.contracts import PabloCrosslinkRequirement
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestor
from polyzymd.builders.conjugation.product_pablo import build_product_state_pablo_library
from polyzymd.config.schema import ConjugationCcdPabloPolicyConfig


def test_product_state_pablo_library_preserves_chain_c_residues(tmp_path: Path):
    """The POC library should define LYX/NHX/SBM without a collapsed POLY residue."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_product_pdb(), encoding="utf-8")
    plan = _resolved_plan_like()

    library = build_product_state_pablo_library(
        product,
        source,
        None,
        _generated_fragment_like(),
        plan,
    )

    summaries = {summary.residue_name: summary for summary in library.summaries}
    assert {"LYX", "NHX", "SBM"}.issubset(summaries)
    assert "POLY" not in summaries
    assert [summary.residue_name for summary in library.summaries if summary.chain_id == "C"] == [
        "NHX",
        "SBM",
        "EGP",
    ]
    assert summaries["LYX"].crosslink == ("NZ", "C047")
    assert summaries["LYX"].linking_bond == ("C", "N")
    assert {"H2", "OXT", "HXT"}.issubset(summaries["LYX"].leaving_atom_names)
    lyx_definition = next(
        definition for definition in library.definitions if definition.residue_name == "LYX"
    )
    lyx_nz = next(atom for atom in lyx_definition.atoms if atom.name == "NZ")
    assert lyx_nz.charge == 0
    assert summaries["NHX"].crosslink == ("C047", "NZ")
    assert {"HZ2", "HZ3"}.issubset(summaries["LYX"].leaving_atom_names)
    assert "LG" in summaries["NHX"].leaving_atom_names
    chain_c = [summary for summary in library.summaries if summary.chain_id == "C"]
    assert [summary.residue_name for summary in chain_c] == ["NHX", "SBM", "EGP"]
    assert chain_c[0].linking_bond == ("POU", "PIN")
    assert chain_c[1].linking_bond == ("POU", "PIN")
    assert chain_c[2].linking_bond == ("POU", "PIN")
    assert chain_c[1].crosslink is None
    assert chain_c[2].crosslink is None
    assert {"L001", "R001", "L002", "R002"}.issubset(
        {atom for summary in chain_c for atom in summary.leaving_atom_names}
    )


def test_ingest_structure_accepts_prebuilt_residue_library(monkeypatch, tmp_path: Path):
    """A caller-supplied library should bypass policy with_crosslink construction."""
    import polyzymd.builders.conjugation.pablo_adapter as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n", encoding="utf-8")
    supplied_library = object()
    received = []

    class FakeTopology:
        n_molecules = 0
        n_bonds = 0

        @property
        def atoms(self):
            return iter(())

        @property
        def bonds(self):
            return iter(())

    def fake_topology_from_pdb(*args, **kwargs):
        received.append(kwargs["residue_library"])
        return FakeTopology()

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=object(),
        topology_from_pdb=fake_topology_from_pdb,
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")
    policy = ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            {
                "residues": ("LYX", "NHX"),
                "linking_atoms": ("NZ", "C047"),
                "leaving_atoms": ((), ()),
            }
        ]
    )

    result = PabloIngestor(policy=policy).ingest_structure(
        structure,
        residue_library=supplied_library,
    )

    assert result.success is True
    assert received == [supplied_library]


def _resolved_plan_like():
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("HZ2", "HZ3"), ("LG",)),
        bond_order=1,
    )
    return SimpleNamespace(
        pablo_crosslink_requirement=requirement,
        protein_link_atom=SimpleNamespace(chain_id="A", residue_number=23),
        modifier_link_atom=SimpleNamespace(chain_id="C", residue_number=5),
        modifier_leaving_atoms=(
            SimpleNamespace(atom_name="LG", element="O", charge="", residue_number=5),
        ),
        contract=SimpleNamespace(
            protein_endpoint=SimpleNamespace(
                selector=SimpleNamespace(residue_name="LYS"),
            ),
        ),
    )


def _generated_fragment_like():
    return SimpleNamespace(
        bonds=(("C047", "O020"), ("C001", "C002")),
        bond_orders=(("C047", "O020", 2), ("C001", "C002", 1)),
    )


def _source_lys_pdb() -> str:
    return "".join(
        [
            _pdb_atom(1, "N", "LYS", "A", 23, "N"),
            _pdb_atom(2, "CA", "LYS", "A", 23, "C"),
            _pdb_atom(3, "C", "LYS", "A", 23, "C"),
            _pdb_atom(4, "O", "LYS", "A", 23, "O"),
            _pdb_atom(5, "CB", "LYS", "A", 23, "C"),
            _pdb_atom(6, "CG", "LYS", "A", 23, "C"),
            _pdb_atom(7, "CD", "LYS", "A", 23, "C"),
            _pdb_atom(8, "CE", "LYS", "A", 23, "C"),
            _pdb_atom(9, "NZ", "LYS", "A", 23, "N"),
            _pdb_atom(10, "HZ1", "LYS", "A", 23, "H"),
            _pdb_atom(11, "HZ2", "LYS", "A", 23, "H"),
            _pdb_atom(12, "HZ3", "LYS", "A", 23, "H"),
            "END\n",
        ]
    )


def _product_pdb() -> str:
    return "".join(
        [
            _pdb_atom(1, "N", "LYX", "A", 23, "N"),
            _pdb_atom(2, "CA", "LYX", "A", 23, "C"),
            _pdb_atom(3, "C", "LYX", "A", 23, "C"),
            _pdb_atom(4, "O", "LYX", "A", 23, "O"),
            _pdb_atom(5, "CB", "LYX", "A", 23, "C"),
            _pdb_atom(6, "CG", "LYX", "A", 23, "C"),
            _pdb_atom(7, "CD", "LYX", "A", 23, "C"),
            _pdb_atom(8, "CE", "LYX", "A", 23, "C"),
            _pdb_atom(9, "NZ", "LYX", "A", 23, "N"),
            _pdb_atom(10, "HZ1", "LYX", "A", 23, "H"),
            _pdb_atom(11, "C047", "NHX", "C", 5, "C", record="HETATM"),
            _pdb_atom(12, "O020", "NHX", "C", 5, "O", record="HETATM"),
            _pdb_atom(13, "C001", "SBM", "C", 6, "C", record="HETATM"),
            _pdb_atom(14, "C002", "SBM", "C", 6, "C", record="HETATM"),
            _pdb_atom(15, "C003", "EGP", "C", 7, "C", record="HETATM"),
            _pdb_atom(16, "C004", "EGP", "C", 7, "C", record="HETATM"),
            "CONECT    9   11\n",
            "CONECT   11    9   12\n",
            "CONECT   12   13\n",
            "CONECT   13   14\n",
            "CONECT   14   15\n",
            "CONECT   15   16\n",
            "END\n",
        ]
    )


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    element: str,
    *,
    record: str = "ATOM",
) -> str:
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id}"
        f"{residue_number:4d}       0.000   0.000   0.000  1.00  0.00          {element:>2}\n"
    )
