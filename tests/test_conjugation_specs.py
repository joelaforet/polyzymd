"""Tests for internal conjugation attachment build specs."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PabloCrosslinkRequirement,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
)
from polyzymd.builders.conjugation._specs import (
    ConjugationFragment,
    attachment_spec_from_generated_polymer_plan,
    attachment_spec_from_moiety_plan,
)
from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


def test_moiety_adapter_preserves_fragment_data_and_generated_adapter(tmp_path: Path):
    """One-residue moieties should retain chemistry metadata in build specs."""
    pdb_path = tmp_path / "nag.pdb"
    sdf_path = tmp_path / "nag.sdf"
    pdb_path.write_text("END\n", encoding="utf-8")
    sdf_path.write_text("", encoding="utf-8")
    moiety = GeneratedMoietyFragment(
        atoms=(
            _atom(0, 1, "C001", "NAG", formal_charge=1),
            _atom(1, 2, "O002", "NAG", element="O", formal_charge=-1),
        ),
        bonds=((1, 2),),
        bond_orders=((1, 2, 1.0),),
        residues=(_residue(0, "NAG", 1),),
        residue_name="NAG",
        name="n_glycan",
        pdb_path=pdb_path,
        sdf_path=sdf_path,
    )
    plan = _resolved_plan(modifier_residue_name="NAG")

    spec = attachment_spec_from_moiety_plan(
        moiety,
        plan,
        attachment_config=SimpleNamespace(name="glycan_1"),
        attachment_index=1,
        reaction_name="n_glycosylation",
    )

    assert spec.attachment_id == "glycan_1"
    assert spec.fragment.source_kind == "moiety"
    assert spec.fragment.atoms == moiety.atoms
    assert spec.fragment.bonds == moiety.bonds
    assert spec.fragment.bond_orders == moiety.bond_orders
    assert spec.fragment.atoms[0].formal_charge == 1
    assert spec.fragment.sidecars == {"pdb": pdb_path, "sdf": sdf_path}
    assert spec.generated_fragment.reactive_atom_serial == 1
    assert spec.generated_fragment.reactive_atom_index == 0
    assert spec.generated_fragment.leaving_atom_serials == (2,)
    assert spec.generated_fragment.leaving_atom_indices == (1,)
    assert spec.source_fragment is moiety


def test_polymer_adapter_preserves_multi_residue_fragment_and_sdf_sidecar(tmp_path: Path):
    """Generated polymer specs should preserve residue count independently of attachments."""
    sdf_path = tmp_path / "polymer.sdf"
    sdf_path.write_text("", encoding="utf-8")
    polymer = GeneratedPolymerFragment(
        atoms=(
            _atom(0, 10, "C1", "SBM", residue_number=10),
            _atom(1, 11, "RC", "NHS", residue_number=11),
            _atom(2, 12, "LG", "NHS", residue_number=11, element="O"),
        ),
        bonds=((10, 11), (11, 12)),
        bond_orders=((10, 11, 1.0), (11, 12, 2.0)),
        residues=(_residue(0, "SBM", 10), _residue(1, "NHS", 11)),
        sequence="AC",
        reactive_atom_serial=11,
        leaving_atom_serials=(12,),
        name="sbm_nhs",
    )
    plan = _resolved_plan(modifier_residue_name="NHS", modifier_atom_name="RC")

    spec = attachment_spec_from_generated_polymer_plan(
        polymer,
        sdf_path,
        plan,
        attachment_config=SimpleNamespace(name="polymer_1"),
        attachment_index=1,
        reaction_name="nhs_lys",
    )

    assert isinstance(spec.fragment, ConjugationFragment)
    assert spec.fragment.source_kind == "polymer"
    assert len(spec.fragment.residues) == 2
    assert spec.fragment.sequence == "AC"
    assert spec.fragment.sidecars == {"sdf": sdf_path}
    assert spec.source_sidecars == {"sdf": sdf_path}
    assert spec.generated_fragment is polymer
    assert spec.fragment.to_generated_polymer_fragment().residues == polymer.residues


def _atom(
    atom_index: int,
    serial: int,
    atom_name: str,
    residue_name: str,
    *,
    residue_number: int = 1,
    element: str = "C",
    formal_charge: int | None = None,
) -> PolymerFragmentAtom:
    return PolymerFragmentAtom(
        atom_index=atom_index,
        serial=serial,
        atom_name=atom_name,
        residue_name=residue_name,
        residue_number=residue_number,
        sequence_index=0,
        x=0.0,
        y=0.0,
        z=0.0,
        element=element,
        formal_charge=formal_charge,
    )


def _residue(sequence_index: int, residue_name: str, residue_number: int) -> PolymerFragmentResidue:
    return PolymerFragmentResidue(
        sequence_index=sequence_index,
        residue_name=residue_name,
        residue_number=residue_number,
    )


def _resolved_plan(
    *,
    modifier_residue_name: str,
    modifier_atom_name: str = "C001",
) -> ResolvedAttachmentPlan:
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="ASN",
        residue_number=42,
        atom_name="ND2",
    )
    modifier_selector = PdbAtomSelector(
        chain_id="C",
        residue_name=modifier_residue_name,
        residue_number=1,
        atom_name=modifier_atom_name,
        atom_serial=1,
        atom_index=0,
    )
    requirement = PabloCrosslinkRequirement(
        residues=("ASX", modifier_residue_name),
        linking_atoms=("ND2", modifier_atom_name),
        leaving_atoms=((), ("O002",)),
        bond_order=1,
    )
    return ResolvedAttachmentPlan(
        contract=ExplicitLinkageContract(
            protein_endpoint=ReactiveEndpoint(
                participant="protein",
                selector=protein_selector,
                product_residue_name="ASX",
            ),
            modifier_endpoint=ReactiveEndpoint(
                participant="modifier",
                selector=modifier_selector,
                product_residue_name=modifier_residue_name,
                leaving_atom_names=("O002",),
            ),
            bond=LinkageBond(
                protein_atom_selector=protein_selector,
                modifier_atom_selector=modifier_selector,
            ),
            mechanism_name="n_glycosylation",
        ),
        protein_link_atom=PdbAtomRecord(
            serial=42,
            atom_index=41,
            atom_name="ND2",
            residue_name="ASN",
            chain_id="A",
            residue_number=42,
            x=0.0,
            y=0.0,
            z=0.0,
        ),
        modifier_link_atom=PdbAtomRecord(
            serial=1,
            atom_index=0,
            atom_name=modifier_atom_name,
            residue_name=modifier_residue_name,
            chain_id="C",
            residue_number=1,
            x=1.0,
            y=0.0,
            z=0.0,
        ),
        modifier_leaving_atoms=(
            PdbAtomRecord(
                serial=2,
                atom_index=1,
                atom_name="O002",
                residue_name=modifier_residue_name,
                chain_id="C",
                residue_number=1,
                x=1.0,
                y=1.0,
                z=0.0,
            ),
        ),
        protein_product_residue_name="ASX",
        modifier_product_residue_name=modifier_residue_name,
        pablo_crosslink_requirement=requirement,
    )
