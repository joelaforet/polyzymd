"""Tests for product-state Pablo residue-library construction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._linkage import PabloCrosslinkRequirement
from polyzymd.builders.conjugation.pablo import product as product_pablo_module
from polyzymd.builders.conjugation.pablo.charge_templates import build_conjugate_charge_templates
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor
from polyzymd.builders.conjugation.pablo.product import (
    build_product_state_pablo_library,
    build_product_state_pablo_library_for_specs,
)
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.polymer.polymerist import generated_fragment_from_polymerist_pdb
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.config.schema import ConjugationCcdPabloPolicyConfig


@dataclass(frozen=True)
class PolymerChemistryFixture:
    name: str
    atoms: tuple[tuple[str, str, int, str, str], ...]
    bonds: tuple[tuple[str, str, int], ...]
    expected_valences: dict[tuple[int, str], int]


@dataclass(frozen=True)
class _FakeAtomDefinition:
    name: str
    symbol: str = "C"
    charge: int = 0
    leaving: bool = False
    synonyms: tuple[str, ...] = ()

    @classmethod
    def with_defaults(
        cls,
        name: str,
        symbol: str,
        *,
        synonyms: tuple[str, ...] = (),
        charge: int = 0,
        leaving: bool = False,
    ) -> _FakeAtomDefinition:
        """Return a Pablo-like atom definition with defaulted metadata."""
        return cls(name=name, symbol=symbol, charge=charge, leaving=leaving, synonyms=synonyms)

    def replace(self, **updates: object) -> _FakeAtomDefinition:
        """Return a copy with selected Pablo-like fields replaced."""
        values = {
            "name": self.name,
            "symbol": self.symbol,
            "charge": self.charge,
            "leaving": self.leaving,
            "synonyms": self.synonyms,
        }
        values.update(updates)
        return type(self)(**values)


@dataclass(frozen=True)
class _FakeBondDefinition:
    atom1: str
    atom2: str
    order: int = 1

    @classmethod
    def with_defaults(cls, atom1: str, atom2: str, *, order: int = 1) -> _FakeBondDefinition:
        """Return a Pablo-like bond definition with defaulted metadata."""
        return cls(atom1=atom1, atom2=atom2, order=order)


@dataclass(frozen=True)
class _FakeResidueDefinition:
    residue_name: str
    atoms: tuple[_FakeAtomDefinition, ...]
    bonds: tuple[_FakeBondDefinition, ...]
    description: str = "Test definition"
    linking_bond: _FakeBondDefinition | None = None
    crosslink: _FakeBondDefinition | None = None
    virtual_sites: tuple[object, ...] = ()

    def replace(self, **updates: object) -> _FakeResidueDefinition:
        """Return a copy with selected Pablo-like fields replaced."""
        values = {
            "residue_name": self.residue_name,
            "atoms": self.atoms,
            "bonds": self.bonds,
            "description": self.description,
            "linking_bond": self.linking_bond,
            "crosslink": self.crosslink,
            "virtual_sites": self.virtual_sites,
        }
        values.update(updates)
        return type(self)(**values)


_SBMA_EGPMA_NHS_CHEMISTRY = PolymerChemistryFixture(
    name="sbma_egpma_nhs_explicit_valence",
    atoms=(
        ("CAA", "NHX", 1, "C", ""),
        ("OAA", "NHX", 1, "O", ""),
        ("CBA", "NHX", 1, "C", ""),
        ("CCA", "SBM", 2, "C", ""),
        ("OCA", "SBM", 2, "O", ""),
        ("OCB", "SBM", 2, "O", ""),
        ("CSB", "SBM", 2, "C", ""),
        ("CNQ", "SBM", 2, "C", ""),
        ("NQA", "SBM", 2, "N", "1+"),
        ("M1", "SBM", 2, "C", ""),
        ("M2", "SBM", 2, "C", ""),
        ("M3", "SBM", 2, "C", ""),
        ("SUL", "SBM", 2, "S", ""),
        ("OS1", "SBM", 2, "O", ""),
        ("OS2", "SBM", 2, "O", ""),
        ("OS3", "SBM", 2, "O", "1-"),
        ("E1", "EGP", 3, "C", ""),
        ("E2", "EGP", 3, "C", ""),
        ("EO1", "EGP", 3, "O", ""),
    ),
    bonds=(
        ("CAA", "OAA", 2),
        ("CAA", "CBA", 1),
        ("CCA", "OCA", 2),
        ("CCA", "OCB", 1),
        ("CCA", "CSB", 1),
        ("CSB", "CNQ", 1),
        ("CNQ", "NQA", 1),
        ("NQA", "M1", 1),
        ("NQA", "M2", 1),
        ("NQA", "M3", 1),
        ("CSB", "SUL", 1),
        ("SUL", "OS1", 2),
        ("SUL", "OS2", 2),
        ("SUL", "OS3", 1),
        ("E1", "E2", 2),
        ("E2", "EO1", 1),
    ),
    expected_valences={
        (1, "CAA"): 4,
        (2, "CCA"): 4,
        (2, "NQA"): 4,
        (2, "SUL"): 6,
        (3, "E2"): 3,
    },
)


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


def test_product_state_guard_allows_single_residue_smiles_moiety():
    """One-residue SMILES moieties such as NAG should not look like polymer collapse."""
    summaries = [
        product_pablo_module.ProductStatePabloDefinitionSummary(
            residue_name="ASX",
            chain_id="B",
            residue_number=625,
        ),
        product_pablo_module.ProductStatePabloDefinitionSummary(
            residue_name="NAG",
            chain_id="C",
            residue_number=1,
        ),
    ]

    product_pablo_module._validate_no_whole_polymer_collapse(summaries)


def test_product_state_pablo_library_uses_connectivity_for_chain_c_links(tmp_path: Path):
    """Polymer links should not require adjacent residue numbers in emitted PDB order."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_public_three_mer_product_pdb(), encoding="utf-8")
    plan = _public_three_mer_resolved_plan_like()

    library = build_product_state_pablo_library(
        product,
        source,
        None,
        None,
        plan,
    )

    chain_c = [summary for summary in library.summaries if summary.chain_id == "C"]
    assert [summary.residue_name for summary in chain_c] == ["NHX", "SBM", "EGP"]
    assert all(summary.linking_bond == ("POU", "PIN") for summary in chain_c)
    assert {"PIN", "POU"}.issubset(chain_c[0].atom_names)
    assert {"PIN", "POU"} & set(chain_c[2].atom_names)
    assert any("chain-C connectivity" in diagnostic for diagnostic in library.diagnostics)


def test_ingest_structure_accepts_prebuilt_residue_library(monkeypatch, tmp_path: Path):
    """A caller-supplied library should bypass policy with_crosslink construction."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

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


def test_product_state_pablo_library_for_specs_uses_generic_fragments_and_sidecars(
    monkeypatch,
    tmp_path: Path,
):
    """Plural product library generation should consume spec sidecars, not moiety fields."""
    product = tmp_path / "product.pdb"
    source = tmp_path / "source.pdb"
    sdf_paths = (tmp_path / "one.sdf", tmp_path / "two.sdf")
    for path in (product, source, *sdf_paths):
        path.write_text("END\n", encoding="utf-8")
    requirement = PabloCrosslinkRequirement(
        residues=("ASX", "NAG"),
        linking_atoms=("ND2", "C001"),
        leaving_atoms=((), ()),
        bond_order=1,
    )
    specs = tuple(
        SimpleNamespace(
            attachment_id=f"spec_{index}",
            fragment=object(),
            source_sidecars={"sdf": sdf_path},
            pablo_crosslink_requirement=requirement,
        )
        for index, sdf_path in enumerate(sdf_paths, start=1)
    )
    calls = []

    def fake_singular(**kwargs):
        calls.append(kwargs)
        return product_pablo_module.ProductStatePabloLibrary(
            residue_library=object(),
            definitions=(SimpleNamespace(name=f"definition_{len(calls)}"),),
            summaries=(),
            crosslink_requirement=requirement,
        )

    fake_pablo = SimpleNamespace(
        STD_CCD_CACHE=SimpleNamespace(with_=lambda definitions: ("combined", definitions))
    )
    monkeypatch.setattr(product_pablo_module, "build_product_state_pablo_library", fake_singular)

    library = build_product_state_pablo_library_for_specs(
        product,
        source,
        specs,
        pablo_module=fake_pablo,
    )

    assert [call["polymer_sdf"] for call in calls] == list(sdf_paths)
    assert [call["generated_fragment"] for call in calls] == [spec.fragment for spec in specs]
    assert [call["resolved_plan"] for call in calls] == list(specs)
    assert library.residue_library[0] == "combined"
    assert len(library.definitions) == 2


@pytest.mark.parametrize(
    "fixture",
    [_SBMA_EGPMA_NHS_CHEMISTRY],
    ids=lambda fixture: fixture.name,
)
def test_product_state_pablo_library_preserves_fixture_bond_orders_and_valence(
    fixture: PolymerChemistryFixture,
    tmp_path: Path,
):
    """Product-state definitions should preserve explicit input chemistry generically."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture), encoding="utf-8")
    polymer_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")

    library = build_product_state_pablo_library(
        product,
        source,
        polymer_sdf,
        _generated_fragment_from_fixture(fixture),
        _fixture_resolved_plan_like(),
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert [summary.residue_name for summary, _definition in _chain_c_definitions(library)] == [
        "NHX",
        "SBM",
        "EGP",
    ]
    for residue_number, atom_name in fixture.expected_valences:
        assert (
            _definition_valence(definitions[residue_number], atom_name)
            == fixture.expected_valences[(residue_number, atom_name)]
        )
    for atom1, atom2, order in fixture.bonds:
        residue_number = _fixture_residue_number(fixture, atom1)
        assert _bond_order(definitions[residue_number], atom1, atom2) == order
    assert _definition_atom_charge(definitions[2], "NQA") == 1
    assert _definition_atom_charge(definitions[2], "OS3") == -1
    assert library.residue_partial_charges == ()
    with pytest.raises(ValueError, match="no validated/explicit partial-charge provenance"):
        build_conjugate_charge_templates(SimpleNamespace(molecules=()), library)


@pytest.mark.parametrize(
    "fixture",
    [_SBMA_EGPMA_NHS_CHEMISTRY],
    ids=lambda fixture: fixture.name,
)
def test_product_state_pablo_library_maps_sdf_orders_after_product_residue_renumbering(
    fixture: PolymerChemistryFixture,
    tmp_path: Path,
):
    """SDF bond orders should map by unique atom identity after product residue renumbering."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(
        _fixture_product_pdb(fixture, residue_number_map={1: 3, 2: 1, 3: 2}),
        encoding="utf-8",
    )
    polymer_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = _generated_fragment_from_fixture(fixture, include_bond_orders=False)

    library = build_product_state_pablo_library(
        product,
        source,
        polymer_sdf,
        fragment,
        _fixture_resolved_plan_like(modifier_residue_number=3),
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert _bond_order(definitions[1], "SUL", "OS1") == 2
    assert _bond_order(definitions[1], "SUL", "OS2") == 2
    assert _bond_order(definitions[3], "CAA", "OAA") == 2


def test_product_state_pablo_library_preserves_sdf_formal_charges_without_pdb_charges(
    tmp_path: Path,
):
    """Product definitions should use generated-fragment SDF charges before PDB fallback."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_pdb = tmp_path / "polymer.pdb"
    polymer_sdf = polymer_pdb.with_suffix(".sdf")
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture, include_charges=False), encoding="utf-8")
    polymer_pdb.write_text(_fixture_polymer_pdb(fixture, include_charges=False), encoding="utf-8")
    polymer_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = generated_fragment_from_polymerist_pdb(
        polymer_pdb,
        reactive_atom_name="CAA",
        leaving_atom_names=("OAA",),
    )

    library = build_product_state_pablo_library(
        product,
        source,
        polymer_sdf,
        fragment,
        _fixture_resolved_plan_like(),
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert _definition_atom_charge(definitions[2], "NQA") == 1
    assert _definition_atom_charge(definitions[2], "OS3") == -1


def test_product_state_pablo_library_uses_charged_sdf_formal_charges(tmp_path: Path):
    """Product formal charges should come from charged SDF sidecars when present."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    neutral_atoms = tuple(
        (atom_name, residue_name, residue_number, element, "")
        for atom_name, residue_name, residue_number, element, _charge in fixture.atoms
    )
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_pdb = tmp_path / "polymer.pdb"
    raw_sdf = tmp_path / "polymer_raw.sdf"
    charged_sdf = tmp_path / "polymer_charged.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture, include_charges=False), encoding="utf-8")
    polymer_pdb.write_text(_fixture_polymer_pdb(fixture, include_charges=False), encoding="utf-8")
    polymer_pdb.with_suffix(".sdf").write_text(
        _sdf_from_atoms_and_bonds(neutral_atoms, fixture.bonds),
        encoding="utf-8",
    )
    raw_sdf.write_text(_sdf_from_atoms_and_bonds(neutral_atoms, fixture.bonds), encoding="utf-8")
    charged_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = generated_fragment_from_polymerist_pdb(
        polymer_pdb,
        reactive_atom_name="CAA",
        leaving_atom_names=("OAA",),
    )

    library = build_product_state_pablo_library(
        product,
        source,
        raw_sdf,
        fragment,
        _fixture_resolved_plan_like(),
        charged_polymer_sdf=charged_sdf,
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert _definition_atom_charge(definitions[2], "NQA") == 1
    assert _definition_atom_charge(definitions[2], "OS3") == -1


def test_product_state_pablo_library_maps_charged_sdf_by_atom_index(tmp_path: Path):
    """Charged SDF formal charges should preserve atom-index mapping."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    neutral_atoms = tuple(
        (atom_name, residue_name, residue_number, element, "")
        for atom_name, residue_name, residue_number, element, _charge in fixture.atoms
    )
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    raw_sdf = tmp_path / "polymer_raw.sdf"
    charged_sdf = tmp_path / "polymer_charged.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture, include_charges=False), encoding="utf-8")
    raw_sdf.write_text(_sdf_from_atoms_and_bonds(neutral_atoms, fixture.bonds), encoding="utf-8")
    charged_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = _generated_fragment_from_fixture(fixture, include_bond_orders=False)
    fragment = fragment.model_copy(update={"atoms": tuple(reversed(fragment.atoms))})

    library = build_product_state_pablo_library(
        product,
        source,
        raw_sdf,
        fragment,
        _fixture_resolved_plan_like(),
        charged_polymer_sdf=charged_sdf,
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert _definition_atom_charge(definitions[2], "NQA") == 1
    assert _definition_atom_charge(definitions[2], "OS3") == -1


def test_product_state_pablo_library_for_specs_scopes_repeated_polymer_residues(
    tmp_path: Path,
):
    """Repeated NHS-polymer fragments should map bond orders to each emitted copy."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    source.write_text(_source_two_lys_pdb(), encoding="utf-8")
    product.write_text(_two_attachment_fixture_product_pdb(fixture), encoding="utf-8")

    fragment = _generated_fragment_from_fixture(fixture)
    specs = (
        SimpleNamespace(
            attachment_id="site_23",
            fragment=fragment,
            **_fixture_resolved_plan_like(
                protein_residue_number=23, modifier_residue_number=6
            ).__dict__,
            product_residue_mappings=_product_residue_mappings((1, 2, 3), (6, 7, 8)),
        ),
        SimpleNamespace(
            attachment_id="site_44",
            fragment=fragment,
            **_fixture_resolved_plan_like(
                protein_residue_number=44, modifier_residue_number=16
            ).__dict__,
            product_residue_mappings=_product_residue_mappings((1, 2, 3), (16, 17, 18)),
        ),
    )

    library = build_product_state_pablo_library_for_specs(product, source, specs)

    chain_c = _chain_c_definitions(library)
    assert [summary.residue_number for summary, _definition in chain_c] == [6, 7, 8, 16, 17, 18]
    definitions = {summary.residue_number: definition for summary, definition in chain_c}
    for residue_number in (7, 17):
        assert _definition_valence(definitions[residue_number], "SUL") == 6
        assert _definition_atom_charge(definitions[residue_number], "NQA") == 1
        assert _definition_atom_charge(definitions[residue_number], "OS3") == -1
    for residue_number in (6, 16):
        assert _bond_order(definitions[residue_number], "CAA", "OAA") == 2
        assert _bond_order(definitions[residue_number], "CAA", "CBA") == 1


def test_product_state_pablo_library_maps_repeated_residue_name_bond_orders(
    tmp_path: Path,
):
    """Repeated canonical glycan residues should receive source-local bond orders."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_repeated_glycan_product_pdb(), encoding="utf-8")
    fragment = _repeated_glycan_fragment()
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "4YB"),
        linking_atoms=("NZ", "C1"),
        leaving_atoms=(("HZ2", "HZ3"), ()),
        bond_order=1,
    )
    resolved_plan = SimpleNamespace(
        pablo_crosslink_requirement=requirement,
        protein_link_atom=SimpleNamespace(chain_id="A", residue_number=23),
        modifier_link_atom=SimpleNamespace(chain_id="C", residue_number=10),
        modifier_leaving_atoms=(),
        contract=SimpleNamespace(
            protein_endpoint=SimpleNamespace(
                selector=SimpleNamespace(residue_name="LYS"),
            ),
        ),
    )

    library = build_product_state_pablo_library(
        product,
        source,
        None,
        fragment,
        resolved_plan,
        product_residue_mappings=_product_residue_mappings((2, 3), (10, 11)),
    )

    definitions = {
        summary.residue_number: definition for summary, definition in _chain_c_definitions(library)
    }
    assert _bond_order(definitions[10], "C2N", "O2N") == 2
    assert _bond_order(definitions[11], "C2N", "O2N") == 2


def test_product_state_pablo_library_for_specs_deduplicates_five_repeated_templates(
    tmp_path: Path,
):
    """Repeated generated crosslink labels should share one internal Pablo template."""
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    raw_crosslink_atoms = ("C063", "C071", "C079", "C087", "C095")
    residue_numbers = (23, 44, 65, 86, 107)
    source.write_text(_source_many_lys_pdb(residue_numbers), encoding="utf-8")
    product.write_text(
        _repeated_generated_attachment_product_pdb(residue_numbers, raw_crosslink_atoms),
        encoding="utf-8",
    )
    specs = tuple(
        SimpleNamespace(
            attachment_id=f"site_{residue_number}",
            fragment=SimpleNamespace(bonds=((raw_atom, "O020"),), bond_orders=()),
            **_generated_crosslink_resolved_plan_like(
                protein_residue_number=residue_number,
                modifier_residue_number=index,
                modifier_link_atom=raw_atom,
            ).__dict__,
            product_residue_mappings=_product_residue_mappings((1,), (index,)),
        )
        for index, (residue_number, raw_atom) in enumerate(
            zip(residue_numbers, raw_crosslink_atoms, strict=True), start=1
        )
    )

    library = build_product_state_pablo_library_for_specs(product, source, specs)

    assert [definition.residue_name for definition in library.definitions].count("LYX") == 1
    assert [definition.residue_name for definition in library.definitions].count("NHX") == 1
    nhx_definition = next(
        definition for definition in library.definitions if definition.residue_name == "NHX"
    )
    assert _definition_atom_synonyms(nhx_definition, "CXL") == raw_crosslink_atoms
    assert _definition_crosslink(nhx_definition) == ("CXL", "NZ")
    lyx_definition = next(
        definition for definition in library.definitions if definition.residue_name == "LYX"
    )
    assert _definition_crosslink(lyx_definition) == ("CXL", "NZ")
    assert all("CXL" not in summary.atom_names for summary in library.summaries)
    assert all(
        summary.crosslink != ("CXL", "NZ") and summary.crosslink != ("NZ", "CXL")
        for summary in library.summaries
    )
    assert library.crosslink_requirement.linking_atoms == ("NZ", "C063")


def test_product_state_pablo_definitions_allow_sb2_variants_with_distinct_selectors():
    """Same residue names may coexist when non-leaving atom selectors are unique."""
    first = _fake_residue_definition("SB2", ("C1", "O1"), (("C1", "O1", 1),))
    second = _fake_residue_definition("SB2", ("C2", "O2"), (("C2", "O2", 1),))

    assert product_pablo_module._deduplicate_product_definitions((first, second)) == (first, second)


def test_product_state_pablo_definitions_allow_distinct_leaving_selectors():
    """Same retained atoms may coexist when leaving selectors are distinct."""
    first = _fake_residue_definition(
        "SB2",
        ("C1", "O1"),
        (("C1", "O1", 1), ("C1", "LG1", 1)),
        leaving_atom_names=("LG1",),
    )
    second = _fake_residue_definition(
        "SB2",
        ("C1", "O1"),
        (("C1", "O1", 1), ("C1", "LG2", 1)),
        leaving_atom_names=("LG2",),
    )

    assert product_pablo_module._deduplicate_product_definitions((first, second)) == (
        first,
        second,
    )


def test_product_state_pablo_coalesces_linking_atom_leaving_fragments():
    """Polymer link definitions should expose one direct leaving fragment per link atom."""
    definition = product_pablo_module._build_pdb_residue_definition(
        _FakeAtomDefinition,
        _FakeBondDefinition,
        _FakeResidueDefinition,
        residue_atoms=(
            _pdb_record(1, "POU", "VMA", residue_number=1),
            _pdb_record(2, "C1", "VMA", residue_number=1),
        ),
        bonds=(("POU", "C1", 1),),
        formal_charges={},
        linking_bond=("POU", "PIN", 1),
        leaving_atoms=(),
        crosslink=None,
        extra_leaving_bonds=(("POU", "L001"), ("POU", "L002")),
        atom_name_aliases={},
        bond_order=1,
        description="Test product-state definition",
    )

    leaving_names = {atom.name for atom in definition.atoms if atom.leaving}
    direct_leavers = {
        bond.atom2 if bond.atom1 == "POU" else bond.atom1
        for bond in definition.bonds
        if "POU" in (bond.atom1, bond.atom2)
        and (bond.atom2 if bond.atom1 == "POU" else bond.atom1) in leaving_names
    }

    assert direct_leavers == {"L001"}
    assert _FakeBondDefinition("L001", "L002", 1) in definition.bonds


def test_product_state_pablo_plans_degree_three_graph_with_one_crosslink():
    """Degree-three glycan nodes should use prior/posterior plus one Pablo crosslink."""
    root = ("C", "4YB", 1, "")
    middle = ("C", "0YB", 2, "")
    posterior = ("C", "0MB", 3, "")
    branch = ("C", "1MA", 4, "")
    plans, diagnostics = product_pablo_module._plan_polymer_external_links(
        (
            product_pablo_module._PolymerExternalBond(root, middle, "O4", "C1", 1),
            product_pablo_module._PolymerExternalBond(middle, posterior, "O4", "C1", 1),
            product_pablo_module._PolymerExternalBond(middle, branch, "O3", "C1", 1),
        ),
        reserved_crosslink_keys={root},
    )

    assert any("Branched product polymer connectivity" in item for item in diagnostics)
    assert plans[root].linking_bond == ("POU", "PIN", 1)
    assert plans[posterior].linking_bond == ("POU", "PIN", 1)
    assert plans[middle].crosslink == ("O3", "C1")
    assert plans[branch].crosslink == ("C1", "O3")


def test_product_state_pablo_keeps_pendant_branch_subtrees_linked():
    """Pendant branch subtrees should not leave child edges as unmatched CONECT records."""
    keys = {index: ("C", f"Z{index:02d}", index, "") for index in range(1, 9)}

    plans, _diagnostics = product_pablo_module._plan_polymer_external_links(
        (
            product_pablo_module._PolymerExternalBond(keys[1], keys[2], "O4", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[2], keys[3], "O4", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[3], keys[4], "O6", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[4], keys[5], "O6", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[4], keys[6], "O3", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[3], keys[7], "O3", "C1", 1),
            product_pablo_module._PolymerExternalBond(keys[7], keys[8], "O2", "C1", 1),
        ),
        reserved_crosslink_keys={keys[1]},
    )

    assert plans[keys[7]].linking_bond == ("POU", "PIN", 1)
    assert plans[keys[8]].linking_bond == ("POU", "PIN", 1)
    assert plans[keys[7]].crosslink == ("C1", "O3")
    assert plans[keys[8]].crosslink is None


def test_product_state_pablo_prefers_conect_external_edges_over_fragment_aliases():
    """PDB CONECT should be authoritative for scoped inter-residue glycan edges."""
    product_atoms = [
        _pdb_record(1, "O4", "Z01", residue_number=1),
        _pdb_record(2, "O3", "Z01", residue_number=1),
        _pdb_record(3, "C1", "Z02", residue_number=2),
        _pdb_record(4, "H4O", "Z01", residue_number=1),
    ]

    _product_bonds, external_bonds = product_pablo_module._product_bonds_and_links_by_residue(
        product_atoms,
        ((1, 3),),
        fragment_bonds=(
            ((1, "O3"), (2, "C1"), 1),
            ((1, "H4O"), (2, "C1"), 1),
        ),
    )

    assert external_bonds == (
        product_pablo_module._PolymerExternalBond(
            ("C", "Z01", 1, ""),
            ("C", "Z02", 2, ""),
            "O4",
            "C1",
            1,
        ),
    )


def test_product_state_pablo_restores_glycam_n_acetyl_carbonyl_order():
    """GLYCAM N-acetyl carbonyls should not become radical single-bond motifs."""
    product_atoms = [
        _pdb_record(1, "C2N", "Z04", residue_number=4),
        _pdb_record(2, "O2N", "Z04", residue_number=4),
        _pdb_record(3, "N2", "Z04", residue_number=4),
        _pdb_record(4, "CME", "Z04", residue_number=4),
    ]

    product_bonds, _external_bonds = product_pablo_module._product_bonds_and_links_by_residue(
        product_atoms,
        ((1, 2), (1, 3), (1, 4)),
        fragment_bonds=(),
    )

    assert product_bonds[("C", "Z04", 4, "")][0] == ("C2N", "CME", 1)
    assert ("C2N", "O2N", 2) in product_bonds[("C", "Z04", 4, "")]


def test_product_state_pablo_legacy_direct_serial_fallback_requires_unambiguous_name():
    """Legacy direct specs may fallback from source serials only when name matching is unique."""
    diagnostics: list[str] = []
    atoms = [
        _pdb_record(4, "CB", "ASX", residue_number=1),
        _pdb_record(10, "C1", "Z01", residue_number=2),
    ]

    key = product_pablo_module._locate_residue_key(
        atoms,
        residue_name="Z01",
        atom_name="C1",
        resolved_atom=SimpleNamespace(serial=4),
        diagnostics=diagnostics,
        allow_legacy_serial_fallback=True,
    )

    assert key == ("C", "Z01", 2, "")
    assert any(
        "Legacy direct product-state endpoint serial fallback" in item for item in diagnostics
    )


def test_product_state_pablo_rejects_degree_greater_than_three_graph():
    """Degree-four glycan graphs should fail before Pablo emits opaque errors."""
    hub = ("C", "0YB", 2, "")
    edges = tuple(
        product_pablo_module._PolymerExternalBond(
            hub,
            ("C", f"B{index}", index, ""),
            f"O{index}",
            "C1",
            1,
        )
        for index in range(1, 5)
    )

    with pytest.raises(ValueError, match="degree 4"):
        product_pablo_module._plan_polymer_external_links(
            edges,
            reserved_crosslink_keys={("C", "4YB", 1, "")},
        )


def test_product_state_pablo_rejects_branch_on_reserved_crosslink_residue():
    """A reducing residue cannot spend a second Pablo crosslink on a branch."""
    root = ("C", "4YB", 1, "")
    first = ("C", "0YB", 2, "")
    branch = ("C", "1MA", 3, "")
    second_branch = ("C", "2MA", 4, "")

    with pytest.raises(ValueError, match="reserved protein crosslink"):
        product_pablo_module._plan_polymer_external_links(
            (
                product_pablo_module._PolymerExternalBond(root, first, "O4", "C1", 1),
                product_pablo_module._PolymerExternalBond(root, branch, "O3", "C1", 1),
                product_pablo_module._PolymerExternalBond(root, second_branch, "O6", "C1", 1),
            ),
            reserved_crosslink_keys={root},
        )


def test_product_state_pablo_crosslink_alias_rejects_existing_cxl_name():
    """Generated C### aliases should not mask a real CXL atom in the product residue."""
    residue_atoms = (
        _pdb_record(11, "C063", "NHX", residue_number=1),
        _pdb_record(12, "CXL", "NHX", residue_number=1),
    )

    with pytest.raises(ValueError, match="already contains a real atom"):
        product_pablo_module._modifier_crosslink_atom_aliases(
            "C063",
            residue_atoms=residue_atoms,
            existing_aliases={},
        )


def test_product_state_pablo_crosslink_alias_rejects_duplicate_raw_atom():
    """Generated C### aliases should require one residue-scoped raw PDB atom."""
    residue_atoms = (
        _pdb_record(11, "C063", "NHX", residue_number=1),
        _pdb_record(12, "C063", "NHX", residue_number=1),
    )

    with pytest.raises(ValueError, match="exactly one raw product atom"):
        product_pablo_module._modifier_crosslink_atom_aliases(
            "C063",
            residue_atoms=residue_atoms,
            existing_aliases={},
        )


def test_product_state_pablo_library_rejects_under_specified_sdf(tmp_path: Path):
    """SDF zero-order bonds should fail with an actionable error, not chemistry repair."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture), encoding="utf-8")
    polymer_sdf.write_text(_fixture_sdf(fixture, overrides={("CAA", "OAA"): 0}), encoding="utf-8")

    with pytest.raises(ValueError, match="under-specified zero/unknown bond orders"):
        build_product_state_pablo_library(
            product,
            source,
            polymer_sdf,
            _generated_fragment_from_fixture(fixture),
            _fixture_resolved_plan_like(),
        )


def test_product_state_pablo_library_rejects_sdf_element_order_mismatch(tmp_path: Path):
    """SDF atom order should be validated before bond-order transfer."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    mismatched_atoms = ((fixture.atoms[0][0], *fixture.atoms[0][1:3], "O", ""), *fixture.atoms[1:])
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture), encoding="utf-8")
    polymer_sdf.write_text(
        _sdf_from_atoms_and_bonds(mismatched_atoms, fixture.bonds), encoding="utf-8"
    )

    with pytest.raises(ValueError, match="element/order does not match"):
        build_product_state_pablo_library(
            product,
            source,
            polymer_sdf,
            _generated_fragment_from_fixture(fixture),
            _fixture_resolved_plan_like(),
        )


def test_product_state_pablo_library_rejects_partial_fragment_atom_indices(tmp_path: Path):
    """SDF mapping should require every generated-fragment atom to carry atom_index."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture), encoding="utf-8")
    polymer_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = _generated_fragment_from_fixture(fixture)
    atoms = list(fragment.atoms)
    atoms[3] = atoms[3].model_copy(update={"atom_index": None})
    fragment = fragment.model_copy(update={"atoms": tuple(atoms)})

    with pytest.raises(ValueError, match="complete atom_index values"):
        build_product_state_pablo_library(
            product,
            source,
            polymer_sdf,
            fragment,
            _fixture_resolved_plan_like(),
        )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({0: 1}, "Duplicate atom_index"),
        ({0: -1}, "non-negative integer atom_index"),
        ({0: len(_SBMA_EGPMA_NHS_CHEMISTRY.atoms)}, "out-of-range"),
        ({0: len(_SBMA_EGPMA_NHS_CHEMISTRY.atoms), -1: 0}, "contiguous and exactly match"),
    ],
    ids=("duplicate", "negative", "out_of_range", "non_contiguous"),
)
def test_product_state_pablo_library_rejects_invalid_sdf_bond_atom_indices(
    updates: dict[int, int],
    message: str,
    tmp_path: Path,
):
    """Bond-order SDF mapping should reject invalid generated-fragment atom indices."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    polymer_sdf = tmp_path / "polymer.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture), encoding="utf-8")
    polymer_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = _fragment_with_atom_index_updates(
        _generated_fragment_from_fixture(fixture),
        updates,
    )

    with pytest.raises(ValueError, match=message):
        build_product_state_pablo_library(
            product,
            source,
            polymer_sdf,
            fragment,
            _fixture_resolved_plan_like(),
        )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({0: 1}, "Duplicate atom_index"),
        ({0: -1}, "non-negative integer atom_index"),
        ({0: len(_SBMA_EGPMA_NHS_CHEMISTRY.atoms)}, "out-of-range"),
        ({0: len(_SBMA_EGPMA_NHS_CHEMISTRY.atoms), -1: 0}, "contiguous and exactly match"),
    ],
    ids=("duplicate", "negative", "out_of_range", "non_contiguous"),
)
def test_product_state_pablo_library_rejects_invalid_charged_sdf_atom_indices(
    updates: dict[int, int],
    message: str,
    tmp_path: Path,
):
    """Charged/formal-charge SDF mapping should reject invalid atom indices."""
    fixture = _SBMA_EGPMA_NHS_CHEMISTRY
    neutral_atoms = tuple(
        (atom_name, residue_name, residue_number, element, "")
        for atom_name, residue_name, residue_number, element, _charge in fixture.atoms
    )
    source = tmp_path / "source.pdb"
    product = tmp_path / "product.pdb"
    raw_sdf = tmp_path / "polymer_raw.sdf"
    charged_sdf = tmp_path / "polymer_charged.sdf"
    source.write_text(_source_lys_pdb(), encoding="utf-8")
    product.write_text(_fixture_product_pdb(fixture, include_charges=False), encoding="utf-8")
    raw_sdf.write_text(_sdf_from_atoms_and_bonds(neutral_atoms, fixture.bonds), encoding="utf-8")
    charged_sdf.write_text(_fixture_sdf(fixture), encoding="utf-8")
    fragment = _fragment_with_atom_index_updates(
        _generated_fragment_from_fixture(fixture, include_bond_orders=False),
        updates,
    )

    with pytest.raises(ValueError, match=message):
        build_product_state_pablo_library(
            product,
            source,
            raw_sdf,
            fragment,
            _fixture_resolved_plan_like(),
            charged_polymer_sdf=charged_sdf,
        )


def test_generated_fragment_from_polymerist_pdb_preserves_sdf_bond_orders(tmp_path: Path):
    """GeneratedPolymerFragment should carry explicit SDF orders from the sidecar."""
    pdb_path = tmp_path / "modifier.pdb"
    sdf_path = pdb_path.with_suffix(".sdf")
    atoms = (
        ("C1", "MOL", 1, "C", ""),
        ("O1", "MOL", 1, "O", ""),
        ("C2", "MOL", 1, "C", ""),
    )
    pdb_path.write_text(
        "".join(
            _pdb_atom(index, atom_name, residue_name, "C", residue_number, element, record="HETATM")
            for index, (atom_name, residue_name, residue_number, element, _charge) in enumerate(
                atoms, start=1
            )
        )
        + "CONECT    1    2    3\nCONECT    2    1\nCONECT    3    1\nEND\n",
        encoding="utf-8",
    )
    sdf_path.write_text(
        _sdf_from_atoms_and_bonds(
            atoms,
            (("C1", "O1", 2), ("C1", "C2", 1)),
        ),
        encoding="utf-8",
    )

    fragment = generated_fragment_from_polymerist_pdb(
        pdb_path,
        reactive_atom_name="C1",
        leaving_atom_names=("O1",),
    )

    assert _fragment_bond_order(fragment, 1, 2) == 2
    assert _fragment_bond_order(fragment, 1, 3) == 1


def _fixture_resolved_plan_like(
    *,
    protein_residue_number: int = 23,
    modifier_residue_number: int = 1,
):
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "CAA"),
        leaving_atoms=(("HZ2", "HZ3"), ()),
        bond_order=1,
    )
    return SimpleNamespace(
        pablo_crosslink_requirement=requirement,
        protein_link_atom=SimpleNamespace(chain_id="A", residue_number=protein_residue_number),
        modifier_link_atom=SimpleNamespace(chain_id="C", residue_number=modifier_residue_number),
        modifier_leaving_atoms=(),
        contract=SimpleNamespace(
            protein_endpoint=SimpleNamespace(
                selector=SimpleNamespace(residue_name="LYS"),
            ),
        ),
    )


def _generated_crosslink_resolved_plan_like(
    *,
    protein_residue_number: int,
    modifier_residue_number: int,
    modifier_link_atom: str,
):
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", modifier_link_atom),
        leaving_atoms=(("HZ2", "HZ3"), ()),
        bond_order=1,
    )
    return SimpleNamespace(
        pablo_crosslink_requirement=requirement,
        protein_link_atom=SimpleNamespace(chain_id="A", residue_number=protein_residue_number),
        modifier_link_atom=SimpleNamespace(chain_id="C", residue_number=modifier_residue_number),
        modifier_leaving_atoms=(),
        contract=SimpleNamespace(
            protein_endpoint=SimpleNamespace(
                selector=SimpleNamespace(residue_name="LYS"),
            ),
        ),
    )


def _generated_fragment_from_fixture(
    fixture: PolymerChemistryFixture,
    *,
    include_bond_orders: bool = True,
) -> GeneratedPolymerFragment:
    atoms = [
        PdbAtomRecord(
            serial=index,
            atom_index=index - 1,
            atom_name=atom_name,
            residue_name=residue_name,
            chain_id="C",
            residue_number=residue_number,
            x=0.0,
            y=0.0,
            z=0.0,
            element=element,
            charge=charge,
            record_name="HETATM",
        )
        for index, (atom_name, residue_name, residue_number, element, charge) in enumerate(
            fixture.atoms,
            start=1,
        )
    ]
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=tuple((atom1, atom2) for atom1, atom2, _order in fixture.bonds),
        bond_orders=fixture.bonds if include_bond_orders else (),
        reactive_atom_name="CAA",
        name=fixture.name,
    )


def _fragment_with_atom_index_updates(
    fragment: GeneratedPolymerFragment,
    updates: dict[int, int],
) -> GeneratedPolymerFragment:
    """Return a fragment copy with selected atom_index values changed."""
    atoms = list(fragment.atoms)
    for position, atom_index in updates.items():
        atoms[position] = atoms[position].model_copy(update={"atom_index": atom_index})
    return fragment.model_copy(update={"atoms": tuple(atoms)})


def _fixture_product_pdb(
    fixture: PolymerChemistryFixture,
    *,
    include_charges: bool = True,
    residue_number_map: dict[int, int] | None = None,
) -> str:
    residue_number_map = residue_number_map or {}
    serial_by_atom = {
        atom_name: index + 11 for index, (atom_name, *_rest) in enumerate(fixture.atoms)
    }
    lines = [*_product_protein_atoms()]
    lines.extend(
        _pdb_atom(
            serial_by_atom[atom_name],
            atom_name,
            residue_name,
            "C",
            residue_number_map.get(residue_number, residue_number),
            element,
            record="HETATM",
            charge=charge if include_charges else "",
        )
        for atom_name, residue_name, residue_number, element, charge in fixture.atoms
    )
    lines.append(f"CONECT    9{serial_by_atom['CAA']:5d}\n")
    conect: dict[int, set[int]] = {}
    for atom1, atom2, _order in fixture.bonds:
        serial1 = serial_by_atom[atom1]
        serial2 = serial_by_atom[atom2]
        conect.setdefault(serial1, set()).add(serial2)
        conect.setdefault(serial2, set()).add(serial1)
    for serial in sorted(conect):
        bonded = "".join(f"{target:5d}" for target in sorted(conect[serial]))
        lines.append(f"CONECT{serial:5d}{bonded}\n")
    lines.append("END\n")
    return "".join(lines)


def _two_attachment_fixture_product_pdb(fixture: PolymerChemistryFixture) -> str:
    lines = [*_product_two_lys_atoms()]
    first_serial_by_atom = _append_fixture_product_copy(
        lines,
        fixture,
        serial_offset=200,
        residue_number_map={1: 6, 2: 7, 3: 8},
    )
    second_serial_by_atom = _append_fixture_product_copy(
        lines,
        fixture,
        serial_offset=400,
        residue_number_map={1: 16, 2: 17, 3: 18},
    )
    lines.append(f"CONECT    9{first_serial_by_atom['CAA']:5d}\n")
    lines.append(f"CONECT  109{second_serial_by_atom['CAA']:5d}\n")
    _append_fixture_conect(lines, fixture, first_serial_by_atom)
    _append_fixture_conect(lines, fixture, second_serial_by_atom)
    lines.append("END\n")
    return "".join(lines)


def _append_fixture_product_copy(
    lines: list[str],
    fixture: PolymerChemistryFixture,
    *,
    serial_offset: int,
    residue_number_map: dict[int, int],
) -> dict[str, int]:
    serial_by_atom = {
        atom_name: serial_offset + index
        for index, (atom_name, *_rest) in enumerate(fixture.atoms, start=1)
    }
    lines.extend(
        _pdb_atom(
            serial_by_atom[atom_name],
            atom_name,
            residue_name,
            "C",
            residue_number_map.get(residue_number, residue_number),
            element,
            record="HETATM",
            charge=charge,
        )
        for atom_name, residue_name, residue_number, element, charge in fixture.atoms
    )
    return serial_by_atom


def _append_fixture_conect(
    lines: list[str],
    fixture: PolymerChemistryFixture,
    serial_by_atom: dict[str, int],
) -> None:
    conect: dict[int, set[int]] = {}
    for atom1, atom2, _order in fixture.bonds:
        serial1 = serial_by_atom[atom1]
        serial2 = serial_by_atom[atom2]
        conect.setdefault(serial1, set()).add(serial2)
        conect.setdefault(serial2, set()).add(serial1)
    for serial in sorted(conect):
        bonded = "".join(f"{target:5d}" for target in sorted(conect[serial]))
        lines.append(f"CONECT{serial:5d}{bonded}\n")


def _fixture_polymer_pdb(
    fixture: PolymerChemistryFixture,
    *,
    include_charges: bool = True,
) -> str:
    lines = [
        _pdb_atom(
            index,
            atom_name,
            residue_name,
            "C",
            residue_number,
            element,
            record="HETATM",
            charge=charge if include_charges else "",
        )
        for index, (atom_name, residue_name, residue_number, element, charge) in enumerate(
            fixture.atoms,
            start=1,
        )
    ]
    serial_by_atom = {
        atom_name: index for index, (atom_name, *_rest) in enumerate(fixture.atoms, 1)
    }
    conect: dict[int, set[int]] = {}
    for atom1, atom2, _order in fixture.bonds:
        serial1 = serial_by_atom[atom1]
        serial2 = serial_by_atom[atom2]
        conect.setdefault(serial1, set()).add(serial2)
        conect.setdefault(serial2, set()).add(serial1)
    for serial in sorted(conect):
        bonded = "".join(f"{target:5d}" for target in sorted(conect[serial]))
        lines.append(f"CONECT{serial:5d}{bonded}\n")
    lines.append("END\n")
    return "".join(lines)


def _fixture_sdf(
    fixture: PolymerChemistryFixture,
    *,
    overrides: dict[tuple[str, str], int] | None = None,
) -> str:
    return _sdf_from_atoms_and_bonds(fixture.atoms, fixture.bonds, overrides=overrides)


def _sdf_from_atoms_and_bonds(
    atoms: tuple[tuple[str, str, int, str, str], ...],
    bonds: tuple[tuple[str, str, int], ...],
    *,
    overrides: dict[tuple[str, str], int] | None = None,
) -> str:
    overrides = overrides or {}
    atom_index = {atom_name: index for index, (atom_name, *_rest) in enumerate(atoms, start=1)}
    lines = [
        "\n",
        "  PolyzyMD test fixture\n",
        "\n",
        f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000\n",
    ]
    for index, (_atom_name, _residue_name, _residue_number, element, charge) in enumerate(atoms):
        x = float(index) * 1.5
        charge_code = _mdl_charge_code(charge)
        lines.append(
            f"{x:10.4f}{0.0:10.4f}{0.0:10.4f} {element:<3} 0  {charge_code:d}"
            "  0  0  0  0  0  0  0  0  0  0\n"
        )
    for atom1, atom2, order in bonds:
        order = overrides.get((atom1, atom2), overrides.get((atom2, atom1), order))
        lines.append(f"{atom_index[atom1]:3d}{atom_index[atom2]:3d}{order:3d}  0\n")
    lines.append("M  END\n$$$$\n")
    return "".join(lines)


def _mdl_charge_code(charge: str) -> int:
    charge_codes = {"3+": 1, "2+": 2, "1+": 3, "1-": 5, "2-": 6, "3-": 7}
    return charge_codes.get(charge, 0)


def _fixture_residue_number(fixture: PolymerChemistryFixture, atom_name: str) -> int:
    for candidate_name, _residue_name, residue_number, _element, _charge in fixture.atoms:
        if candidate_name == atom_name:
            return residue_number
    raise AssertionError(f"Fixture atom {atom_name} was not found")


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


def _public_three_mer_resolved_plan_like():
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C003"),
        leaving_atoms=(("HZ2", "HZ3"), ("LG",)),
        bond_order=1,
    )
    return SimpleNamespace(
        pablo_crosslink_requirement=requirement,
        protein_link_atom=SimpleNamespace(chain_id="A", residue_number=23),
        modifier_link_atom=SimpleNamespace(chain_id="C", residue_number=1),
        modifier_leaving_atoms=(
            SimpleNamespace(atom_name="LG", element="O", charge="", residue_number=1),
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


def _source_two_lys_pdb() -> str:
    lines = []
    serial = 1
    for residue_number in (23, 44):
        for atom_name, element in (
            ("N", "N"),
            ("CA", "C"),
            ("C", "C"),
            ("O", "O"),
            ("CB", "C"),
            ("CG", "C"),
            ("CD", "C"),
            ("CE", "C"),
            ("NZ", "N"),
            ("HZ1", "H"),
            ("HZ2", "H"),
            ("HZ3", "H"),
        ):
            lines.append(_pdb_atom(serial, atom_name, "LYS", "A", residue_number, element))
            serial += 1
    lines.append("END\n")
    return "".join(lines)


def _source_many_lys_pdb(residue_numbers: tuple[int, ...]) -> str:
    """Return source-PDB text with one LYS residue for each requested number."""
    lines = []
    serial = 1
    for residue_number in residue_numbers:
        for atom_name, element in (
            ("N", "N"),
            ("CA", "C"),
            ("C", "C"),
            ("O", "O"),
            ("CB", "C"),
            ("CG", "C"),
            ("CD", "C"),
            ("CE", "C"),
            ("NZ", "N"),
            ("HZ1", "H"),
            ("HZ2", "H"),
            ("HZ3", "H"),
        ):
            lines.append(_pdb_atom(serial, atom_name, "LYS", "A", residue_number, element))
            serial += 1
    lines.append("END\n")
    return "".join(lines)


def _product_residue_mappings(
    source_numbers: tuple[int, ...],
    target_numbers: tuple[int, ...],
) -> dict[str, dict[str, int | str]]:
    return {
        str(source_number): {
            "source_residue_number": source_number,
            "target_chain": "C",
            "target_residue_number": target_number,
        }
        for source_number, target_number in zip(source_numbers, target_numbers, strict=True)
    }


def _repeated_glycan_fragment() -> GeneratedPolymerFragment:
    """Return a two-residue fragment with duplicated GLYCAM residue and atom names."""
    atoms = []
    serial = 1
    for residue_number in (2, 3):
        for atom_name, element in (("C1", "C"), ("C2N", "C"), ("O2N", "O")):
            atoms.append(
                PdbAtomRecord(
                    serial=serial,
                    atom_index=serial - 1,
                    atom_name=atom_name,
                    residue_name="4YB",
                    chain_id="C",
                    residue_number=residue_number,
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    element=element,
                    record_name="HETATM",
                )
            )
            serial += 1
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=((2, 3), (5, 6)),
        bond_orders=((2, 3, 2.0), (5, 6, 2.0)),
        reactive_atom_serial=1,
        name="repeated_4yb",
    )


def _repeated_glycan_product_pdb() -> str:
    """Return a product PDB whose emitted residue numbers differ from the source."""
    lines = [*_product_protein_atoms()]
    serial = 200
    serial_by_residue_atom = {}
    residue_names = {10: "Z11", 11: "Z22"}
    for residue_number in (10, 11):
        for atom_name, element in (("C1", "C"), ("C2N", "C"), ("O2N", "O")):
            lines.append(
                _pdb_atom(
                    serial,
                    atom_name,
                    residue_names[residue_number],
                    "C",
                    residue_number,
                    element,
                    record="HETATM",
                )
            )
            serial_by_residue_atom[(residue_number, atom_name)] = serial
            serial += 1
    lines.append(f"CONECT    9{serial_by_residue_atom[(10, 'C1')]:5d}\n")
    for residue_number in (10, 11):
        c2n = serial_by_residue_atom[(residue_number, "C2N")]
        o2n = serial_by_residue_atom[(residue_number, "O2N")]
        lines.append(f"CONECT{c2n:5d}{o2n:5d}\n")
        lines.append(f"CONECT{o2n:5d}{c2n:5d}\n")
    lines.append("END\n")
    return "".join(lines)


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


def _repeated_generated_attachment_product_pdb(
    residue_numbers: tuple[int, ...],
    raw_crosslink_atoms: tuple[str, ...],
) -> str:
    """Return product-PDB text with repeated LYX/NHX generated attachments."""
    lines: list[str] = []
    serial = 1
    nhx_serials: list[tuple[int, int, int]] = []
    for index, (residue_number, raw_atom) in enumerate(
        zip(residue_numbers, raw_crosslink_atoms, strict=True), start=1
    ):
        nz_serial = serial + 8
        for atom_name, element in (
            ("N", "N"),
            ("CA", "C"),
            ("C", "C"),
            ("O", "O"),
            ("CB", "C"),
            ("CG", "C"),
            ("CD", "C"),
            ("CE", "C"),
            ("NZ", "N"),
            ("HZ1", "H"),
        ):
            lines.append(_pdb_atom(serial, atom_name, "LYX", "A", residue_number, element))
            serial += 1
        carbon_serial = serial
        oxygen_serial = serial + 1
        lines.append(_pdb_atom(carbon_serial, raw_atom, "NHX", "C", index, "C", record="HETATM"))
        lines.append(_pdb_atom(oxygen_serial, "O020", "NHX", "C", index, "O", record="HETATM"))
        serial += 2
        nhx_serials.append((nz_serial, carbon_serial, oxygen_serial))
    for nz_serial, carbon_serial, oxygen_serial in nhx_serials:
        lines.append(f"CONECT{nz_serial:5d}{carbon_serial:5d}\n")
        lines.append(f"CONECT{carbon_serial:5d}{nz_serial:5d}{oxygen_serial:5d}\n")
        lines.append(f"CONECT{oxygen_serial:5d}{carbon_serial:5d}\n")
    lines.append("END\n")
    return "".join(lines)


def _public_three_mer_product_pdb() -> str:
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
            _pdb_atom(11, "C000", "NHX", "C", 1, "C", record="HETATM"),
            _pdb_atom(12, "C001", "NHX", "C", 1, "C", record="HETATM"),
            _pdb_atom(13, "C003", "NHX", "C", 1, "C", record="HETATM"),
            _pdb_atom(14, "O000", "NHX", "C", 1, "O", record="HETATM"),
            _pdb_atom(15, "C008", "SBM", "C", 2, "C", record="HETATM"),
            _pdb_atom(16, "C009", "SBM", "C", 2, "C", record="HETATM"),
            _pdb_atom(17, "C019", "EGP", "C", 3, "C", record="HETATM"),
            _pdb_atom(18, "C020", "EGP", "C", 3, "C", record="HETATM"),
            "CONECT    9   13\n",
            "CONECT   11   12   15\n",
            "CONECT   12   11   17\n",
            "CONECT   13    9   14\n",
            "CONECT   14   13\n",
            "CONECT   15   11   16\n",
            "CONECT   16   15\n",
            "CONECT   17   12   18\n",
            "CONECT   18   17\n",
            "END\n",
        ]
    )


def _product_protein_atoms() -> list[str]:
    return [
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
    ]


def _product_two_lys_atoms() -> list[str]:
    lines = _product_protein_atoms()
    lines.extend(
        _pdb_atom(serial, atom_name, "LYX", "A", 44, element)
        for serial, atom_name, element in (
            (101, "N", "N"),
            (102, "CA", "C"),
            (103, "C", "C"),
            (104, "O", "O"),
            (105, "CB", "C"),
            (106, "CG", "C"),
            (107, "CD", "C"),
            (108, "CE", "C"),
            (109, "NZ", "N"),
            (110, "HZ1", "H"),
        )
    )
    return lines


def _chain_c_definitions(library):
    if len(library.summaries) == len(library.definitions):
        return [
            (summary, definition)
            for summary, definition in zip(library.summaries, library.definitions, strict=True)
            if summary.chain_id == "C"
        ]
    return [
        (summary, _definition_for_summary(summary, library.definitions))
        for summary in library.summaries
        if summary.chain_id == "C"
    ]


def _definition_for_summary(summary, definitions):
    """Return the deduplicated definition represented by a product summary."""
    summary_names = set(summary.atom_names)
    for definition in definitions:
        if definition.residue_name != summary.residue_name:
            continue
        definition_names = {_definition_raw_atom_name(atom) for atom in definition.atoms}
        if definition_names == summary_names:
            return definition
    raise AssertionError(f"Definition for summary {summary.residue_name} was not found")


def _definition_raw_atom_name(atom) -> str:
    synonyms = tuple(getattr(atom, "synonyms", ()) or ())
    if getattr(atom, "name", None) == "CXL" and synonyms:
        return synonyms[0]
    return atom.name


def _bond_order(definition, atom1: str, atom2: str) -> int:
    pair = {atom1, atom2}
    for bond in definition.bonds:
        if {bond.atom1, bond.atom2} == pair:
            return bond.order
    raise AssertionError(f"Bond {atom1}-{atom2} was not defined")


def _definition_valence(definition, atom_name: str) -> int:
    valence = 0
    for bond in definition.bonds:
        if atom_name in {bond.atom1, bond.atom2}:
            valence += int(bond.order)
    crosslink = getattr(definition, "crosslink", None)
    if crosslink is not None and atom_name in {crosslink.atom1, crosslink.atom2}:
        valence += int(crosslink.order)
    linking_bond = getattr(definition, "linking_bond", None)
    if linking_bond is not None and atom_name in {linking_bond.atom1, linking_bond.atom2}:
        valence += int(linking_bond.order)
    return valence


def _definition_atom_charge(definition, atom_name: str) -> int:
    for atom in definition.atoms:
        if atom.name == atom_name:
            return atom.charge
    raise AssertionError(f"Atom {atom_name} was not defined")


def _definition_atom_synonyms(definition, atom_name: str) -> tuple[str, ...]:
    for atom in definition.atoms:
        if atom.name == atom_name:
            return tuple(atom.synonyms)
    raise AssertionError(f"Atom {atom_name} was not defined")


def _definition_crosslink(definition) -> tuple[str, str] | None:
    crosslink = getattr(definition, "crosslink", None)
    if crosslink is None:
        return None
    return tuple(sorted((crosslink.atom1, crosslink.atom2)))


def _fake_residue_definition(
    residue_name: str,
    atom_names: tuple[str, ...],
    bonds: tuple[tuple[str, str, int], ...],
    *,
    leaving_atom_names: tuple[str, ...] = (),
) -> _FakeResidueDefinition:
    """Return a small Pablo-like residue definition for deduplication tests."""
    all_atom_names = (*atom_names, *leaving_atom_names)
    return _FakeResidueDefinition(
        residue_name=residue_name,
        atoms=tuple(
            _FakeAtomDefinition(atom_name, leaving=atom_name in leaving_atom_names)
            for atom_name in all_atom_names
        ),
        bonds=tuple(_FakeBondDefinition(atom1, atom2, order) for atom1, atom2, order in bonds),
    )


def _fragment_bond_order(fragment: GeneratedPolymerFragment, serial1: int, serial2: int) -> int:
    pair = {serial1, serial2}
    for left, right, order in fragment.bond_orders:
        if {left, right} == pair:
            return int(order)
    raise AssertionError(f"Fragment bond {serial1}-{serial2} was not defined")


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    element: str,
    *,
    record: str = "ATOM",
    charge: str = "",
) -> str:
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id}"
        f"{residue_number:4d}       0.000   0.000   0.000  1.00  0.00          "
        f"{element:>2}{charge:>2}\n"
    )


def _pdb_record(
    serial: int,
    atom_name: str,
    residue_name: str,
    *,
    residue_number: int,
) -> PdbAtomRecord:
    """Return a minimal parsed PDB atom record for unit-level alias tests."""
    return PdbAtomRecord(
        serial=serial,
        atom_index=serial - 1,
        atom_name=atom_name,
        residue_name=residue_name,
        chain_id="C",
        residue_number=residue_number,
        x=0.0,
        y=0.0,
        z=0.0,
        element="C",
        record_name="HETATM",
    )
