"""Tests for the mBuild polymer generator and the star-sited SMILES it reads."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

import polyzymd.builders.polymer_generator as polymer_generator_module
from polyzymd.builders.polymer import PolymerBuilder, canonical_sequence
from polyzymd.builders.polymer_generator import (
    DYNAMIC_CACHE_SCHEMA_VERSION,
    MBUILD,
    MIN_CONTACT_ANGSTROM,
    CapUnit,
    MonomerUnit,
    PolymerGenerationError,
    PolymerGenerator,
    assemble_chain,
    closest_contact,
    find_pierced_rings,
    mbuild_provenance,
    parse_star_smiles,
)

#: Star-sited repeat units used in PolyzyMD studies.
MMA = "[*]CC([*])(C)C(=O)OC"
HEMA = "[*]CC([*])(C)C(=O)OCCO"
MEMA = "[*]CC([*])(C)C(=O)OCCOC"
SBMA = "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
EGPMA = "[*]CC([*])(C)C(=O)OCCOc1ccccc1"
MONOMERS = [MMA, HEMA, MEMA, SBMA, EGPMA]

#: Star-sited caps.
H_CAP = "[*][H]"
METHYL_CAP = "[*]C"

#: The Paper 1 SBMA/EGMA chemistry: repeat units and the vinyl caps that end a
#: chain grown by radical polymerization (a whole monomer that kept its double
#: bond), each cap in its monomer's residue.
SBM = "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
EGM = "[*]CC([*])(C)C(=O)OCCOCCOCCOCCOC"
SBM_VINYL = "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
EGM_VINYL = "[*]C=C(C)C(=O)OCCOCCOCCOCCOC"

#: Methyl methacrylate with its quaternary carbon's handedness set (tacticity letters).
MMA_R = "[*]C[C@@]([*])(C)C(=O)OC"
MMA_S = "[*]C[C@]([*])(C)C(=O)OC"

SAGE = ["openff-2.3.0.offxml"]

#: The mBuild release pixi.toml installs.
MBUILD_COMMIT = "84e0f9b00484df185737942af1b3e16b6acde3db"


def _units(**smiles: str) -> dict[str, MonomerUnit]:
    """Return MonomerUnits keyed by label from ``label=smiles`` pairs.

    Each monomer is named and residue-named after its label (``A`` -> ``AAA``).
    """
    return {
        label: MonomerUnit(label, label, label * 3, parse_star_smiles(s, 2, f"monomer {label!r}"))
        for label, s in smiles.items()
    }


def _cap(smiles: str):
    """Return the parsed cap of ``smiles``."""
    return parse_star_smiles(smiles, 1, "cap")


def _option(label: str, smiles: str, probability: float = 1.0, **extra: str) -> dict:
    """Return one cap option as the config gives it."""
    return {"label": label, "smiles": smiles, "probability": probability, **extra}


#: One fixed hydrogen cap (label H) for an end.
H_OPTIONS = [_option("H", H_CAP)]


# ---------------------------------------------------------------------------
# parse_star_smiles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("smiles", MONOMERS)
def test_the_study_monomers_are_read(smiles: str) -> None:
    """Each repeat unit has two anchors, each with a site hydrogen in the '*' slot."""
    unit = parse_star_smiles(smiles, 2, "monomer 'A'")

    assert unit.smiles == smiles
    assert len(unit.anchors) == 2 and unit.anchors[0] != unit.anchors[1]
    for anchor, hydrogen in zip(unit.anchors, unit.site_hydrogens):
        atom = unit.mol.GetAtomWithIdx(hydrogen)
        assert atom.GetAtomicNum() == 1
        assert [n.GetIdx() for n in atom.GetNeighbors()] == [anchor]
    assert all(atom.GetAtomicNum() != 0 for atom in unit.mol.GetAtoms())


def test_the_anchors_follow_the_smiles_order() -> None:
    """The first '*' is the head side (the CH2 of a methacrylate), the second the tail side."""
    unit = parse_star_smiles(MMA, 2, "monomer 'A'")
    head, tail = (unit.mol.GetAtomWithIdx(i) for i in unit.anchors)

    assert head.GetTotalNumHs(includeNeighbors=True) == 3  # CH2 plus its site hydrogen
    assert tail.GetTotalNumHs(includeNeighbors=True) == 1  # the site hydrogen alone
    assert unit.mol.GetNumAtoms() == 15 + 2


@pytest.mark.parametrize("smiles", [H_CAP, METHYL_CAP])
def test_the_caps_are_read(smiles: str) -> None:
    unit = parse_star_smiles(smiles, 1, "end_groups.head")
    assert len(unit.anchors) == 1 and len(unit.site_hydrogens) == 1


@pytest.mark.parametrize(
    ("smiles", "n_sites", "message"),
    [
        ("CC(C)C(=O)OC", 2, "has 0 '\\*' sites; a repeat unit needs exactly 2"),
        ("[*]CC(C)C(=O)OC", 2, "has 1 '\\*' sites; a repeat unit needs exactly 2"),
        ("[*]CC([*])(C[*])C(=O)OC", 2, "has 3 '\\*' sites; a repeat unit needs exactly 2"),
        ("C", 1, "has 0 '\\*' sites; a cap needs exactly 1"),
        ("[*]C[*]", 1, "has 2 '\\*' sites; a cap needs exactly 1"),
        ("[*]CC([*]", 2, "RDKit cannot read the SMILES"),
        ("[*]=CC([*])C", 2, "bonded to one atom by a single bond"),
        ("[*]=C", 1, "bonded to one atom by a single bond"),
        ("[*]C([*])C(=O)OC", 2, "are on one atom"),
        ("[*][*]CC", 2, "bonded to each other"),
        ("[*+]CC([*])(C)C(=O)OC", 2, "carries a charge or hydrogens"),
    ],
    ids=[
        "monomer_0",
        "monomer_1",
        "monomer_3",
        "cap_0",
        "cap_2",
        "unreadable",
        "double_bond",
        "cap_double_bond",
        "one_atom",
        "star_star",
        "charged_star",
    ],
)
def test_star_smiles_the_build_cannot_use_are_refused(smiles, n_sites, message) -> None:
    """Every refusal names the unit and says how to fix it."""
    with pytest.raises(ValueError, match=message) as error:
        parse_star_smiles(smiles, n_sites, "monomer 'Q'")
    assert str(error.value).startswith("monomer 'Q'")
    assert "fix:" in str(error.value)


# ---------------------------------------------------------------------------
# assemble_chain
# ---------------------------------------------------------------------------


def test_an_mma_hema_chain_has_the_atoms_of_its_units() -> None:
    """ABABA of MMA (15 atoms) and HEMA (19) with two hydrogen caps has 85 atoms."""
    chain, residues = assemble_chain("ABABA", _units(A=MMA, B=HEMA), _cap(H_CAP), _cap(H_CAP))

    assert chain.GetNumAtoms() == 3 * 15 + 2 * 19 + 2 == 85
    assert len(residues) == 85


def test_an_sbma_egpma_chain_keeps_its_charges_and_bond_orders() -> None:
    """AABAB of SBMA and EGPMA with a methyl head and hydrogen tail has 180 atoms."""
    chain, residues = assemble_chain(
        "AABAB", _units(A=SBMA, B=EGPMA), _cap(METHYL_CAP), _cap(H_CAP)
    )

    assert chain.GetNumAtoms() == 3 * 39 + 2 * 29 + 4 + 1 == 180
    charges = [atom.GetFormalCharge() for atom in chain.GetAtoms()]
    assert sum(charges) == 0
    nitrogen = [a for a in chain.GetAtoms() if a.GetSymbol() == "N" and a.GetFormalCharge() == 1]
    oxide = [a for a in chain.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() == -1]
    assert len(nitrogen) == 3 and len(oxide) == 3
    carbonyl = Chem.MolFromSmarts("[CX3]=O")
    sulfonyl = Chem.MolFromSmarts("S=O")
    assert len(chain.GetSubstructMatches(carbonyl)) == 5
    assert len(chain.GetSubstructMatches(sulfonyl)) == 2 * 3
    assert chain.GetRingInfo().NumRings() == 2
    assert sum(atom.GetIsAromatic() for atom in chain.GetAtoms()) == 12


def test_residues_are_contiguous_with_the_caps_in_the_end_units() -> None:
    """Atoms come residue by residue, numbered 1..n; each cap joins its end unit."""
    units = _units(A=MMA, B=HEMA)
    chain, residues = assemble_chain("ABABA", units, _cap(METHYL_CAP), _cap(METHYL_CAP))

    assert residues == sorted(residues)
    assert sorted(set(residues)) == [1, 2, 3, 4, 5]
    sizes = [residues.count(n) for n in range(1, 6)]
    assert sizes == [15 + 4, 19, 15, 19, 15 + 4]


def test_the_head_cap_bonds_to_the_first_star() -> None:
    """The head cap replaces the first '*' of the first unit, the tail cap the second of the last."""
    units = _units(A=MMA)
    chain, _ = assemble_chain("AA", units, _cap("[*]Cl"), _cap(H_CAP))
    expected = Chem.CanonSmiles("ClCC(C)(CC(C)C(=O)OC)C(=O)OC")
    assert Chem.MolToSmiles(Chem.RemoveHs(chain)) == expected

    chain, _ = assemble_chain("AA", units, _cap(H_CAP), _cap(H_CAP))
    expected = Chem.CanonSmiles("COC(=O)C(C)(C)CC(C)C(=O)OC")
    assert Chem.MolToSmiles(Chem.RemoveHs(chain)) == expected


@pytest.mark.parametrize("sequence", ["A", "AB"])
def test_chains_of_one_and_two_units_assemble(sequence: str) -> None:
    chain, residues = assemble_chain(sequence, _units(A=MMA, B=HEMA), _cap(H_CAP), _cap(H_CAP))
    assert set(residues) == set(range(1, len(sequence) + 1))
    assert Chem.MolToSmiles(chain)  # sanitized


def test_caps_with_their_own_residue_come_first_and_last() -> None:
    """A cap that is its own residue is residue 1 (head) or the last residue (tail)."""
    units = _units(A=MMA, B=HEMA)
    _, residues = assemble_chain(
        "ABA", units, _cap(METHYL_CAP), _cap(METHYL_CAP), head_residue=True, tail_residue=True
    )
    assert residues == sorted(residues)
    sizes = [residues.count(n) for n in range(1, 6)]
    assert sizes == [4, 15, 19, 15, 4]

    _, residues = assemble_chain(
        "ABA", units, _cap(METHYL_CAP), _cap(METHYL_CAP), head_residue=True
    )
    assert [residues.count(n) for n in range(1, 5)] == [4, 15, 19, 15 + 4]


# ---------------------------------------------------------------------------
# Geometry checks
# ---------------------------------------------------------------------------


def _benzene_and_ethane(ethane_x: float) -> tuple[Chem.Mol, np.ndarray]:
    """A flat benzene ring in the xy plane and an ethane C-C bond along z at ``x``."""
    mol = Chem.MolFromSmiles("c1ccccc1.CC")
    angles = np.arange(6) * np.pi / 3
    ring = np.column_stack((1.4 * np.cos(angles), 1.4 * np.sin(angles), np.zeros(6)))
    ethane = np.array([[ethane_x, 0.0, -0.77], [ethane_x, 0.0, 0.77]])
    return mol, np.vstack((ring, ethane))


def test_a_bond_through_a_ring_is_found() -> None:
    mol, xyz = _benzene_and_ethane(0.0)
    [(ring, bond)] = find_pierced_rings(mol, xyz)
    assert sorted(ring) == [0, 1, 2, 3, 4, 5]
    assert bond == (6, 7)


def test_a_bond_beside_a_ring_is_not_found() -> None:
    mol, xyz = _benzene_and_ethane(3.0)
    assert find_pierced_rings(mol, xyz) == []


def test_closest_contact_skips_atoms_two_bonds_apart() -> None:
    """Contacts count only atoms more than two bonds apart (here, ring vs ethane)."""
    mol, xyz = _benzene_and_ethane(0.0)
    distance, pair = closest_contact(mol, xyz)
    assert distance == pytest.approx(np.hypot(1.4, 0.77))
    assert 6 in pair or 7 in pair

    mol, xyz = _benzene_and_ethane(3.0)
    distance, pair = closest_contact(mol, xyz)
    assert distance == pytest.approx(np.hypot(1.6, 0.77))
    assert set(pair) & {0} and set(pair) & {6, 7}


def test_closest_contact_is_none_when_every_pair_is_close_in_bonds() -> None:
    mol = Chem.MolFromSmiles("CCC")
    assert closest_contact(mol, np.eye(3)) is None


# ---------------------------------------------------------------------------
# The generator without mBuild
# ---------------------------------------------------------------------------


def _generator(tmp_path: Path, **keys: object) -> PolymerGenerator:
    """An MMA/HEMA generator with hydrogen caps caching in ``tmp_path``."""
    arguments = {"cache_directory": tmp_path, "force_field": SAGE, **keys}
    return PolymerGenerator.from_smiles(
        [("A", "MMA", "MMA", MMA), ("B", "HEMA", "HEM", HEMA)],
        head=H_OPTIONS,
        tail=H_OPTIONS,
        **arguments,
    )


def test_from_smiles_reads_the_monomers_and_caps(tmp_path: Path) -> None:
    generator = _generator(tmp_path)
    hema = generator.monomers["B"]
    assert (hema.label, hema.name, hema.residue_name) == ("B", "HEMA", "HEM")
    assert hema.unit.smiles == HEMA and hema.unit.anchors == (1, 2)
    assert generator.end_groups["head"]["H"] == CapUnit(
        "H", None, None, generator.end_groups["head"]["H"].unit
    )
    assert generator.end_groups["tail"]["H"].unit.smiles == H_CAP
    assert generator.charger_type == "nagl"
    assert generator.nagl_model == "openff-gnn-am1bcc-1.0.0.pt"
    assert _generator(tmp_path, charger_type="AM1BCC").nagl_model is None


def test_from_smiles_names_a_bad_monomer(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="^monomer 'B'"):
        PolymerGenerator.from_smiles(
            [("A", "MMA", "MMA", MMA), ("B", "X", "XXX", "[*]CC")],
            head=H_OPTIONS,
            tail=H_OPTIONS,
            cache_directory=tmp_path,
            force_field=SAGE,
        )


def test_the_sdf_and_its_sidecar_are_named_after_the_sequence(tmp_path: Path) -> None:
    generator = _generator(tmp_path, type_prefix="MH")
    sdf = generator.sdf_path("HABAH")
    assert sdf == tmp_path / "MH_seq=HABAH_5-mer_charged.sdf"
    assert generator.metadata_path(sdf) == tmp_path / "MH_seq=HABAH_5-mer_charged.sdf.metadata.json"


def test_the_cache_metadata_records_what_the_chain_was_built_from(tmp_path: Path) -> None:
    """The sidecar records the monomers of the repeat units and the chain's own two caps."""
    generator = PolymerGenerator.from_smiles(
        [("A", "MMA", "MMA", MMA), ("B", "HEMA", "HEM", HEMA)],
        head=[_option("H", H_CAP, 0.5), _option("M", METHYL_CAP, 0.5, name="methyl")],
        tail=[_option("A", "[*]C=C(C)C(=O)OC", residue_name="MMV")],
        cache_directory=tmp_path,
        force_field=SAGE,
    )
    metadata = generator.cache_metadata("MAAAA")

    assert metadata["schema_version"] == DYNAMIC_CACHE_SCHEMA_VERSION == 5
    assert metadata["sequence"] == "MAAAA"
    # The tail cap's label A is not a repeat unit, so HEMA is not recorded.
    assert metadata["monomers"] == {"A": {"name": "MMA", "residue_name": "MMA", "smiles": MMA}}
    assert metadata["end_groups"] == {
        "head": {"label": "M", "name": "methyl", "residue_name": None, "smiles": METHYL_CAP},
        "tail": {
            "label": "A",
            "name": None,
            "residue_name": "MMV",
            "smiles": "[*]C=C(C)C(=O)OC",
        },
    }
    assert metadata["relaxation_force_field"] == SAGE
    assert metadata["builder"]["package"] == "mbuild"
    assert json.loads(json.dumps(metadata)) == metadata


def test_a_sequence_with_unknown_labels_is_refused(tmp_path: Path) -> None:
    with pytest.raises(PolymerGenerationError, match="monomers \\['C'\\]"):
        _generator(tmp_path).build_molecule("HACAH")


def test_split_returns_the_caps_and_the_repeat_units(tmp_path: Path) -> None:
    generator = _generator(tmp_path)
    head, middle, tail = generator.split("HABAH")
    assert (head.label, middle, tail.label) == ("H", "ABA", "H")
    assert head.unit.smiles == tail.unit.smiles == H_CAP


@pytest.mark.parametrize(
    ("sequence", "message"),
    [
        ("HA", "needs a head cap, at least one repeat unit and a tail cap"),
        ("", "needs a head cap"),
        ("XAAH", "head cap 'X'"),
        ("HAAX", "tail cap 'X'"),
        ("HACH", "monomers \\['C'\\]"),
        ("AAAA", "head cap 'A', tail cap 'A'"),
    ],
)
def test_split_refuses_short_sequences_and_unknown_labels(
    tmp_path: Path, sequence: str, message: str
) -> None:
    with pytest.raises(PolymerGenerationError, match=message):
        _generator(tmp_path).split(sequence)


def test_a_build_that_fails_the_checks_is_built_again_with_a_new_seed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A pierced ring or a close contact retries with the next seed, up to max_retries."""
    seeds = []

    def fake_chain(sequence, monomers, head, tail, seed):
        seeds.append(seed)
        return None, None, None

    def fake_coordinates(chain, elements, bonds, xyz):
        return np.random.default_rng(0).random((chain.GetNumAtoms(), 3)) * 10

    class Charger:
        def charge_molecule(self, molecule):
            return molecule

    checks = iter([[((0, 1, 2), (3, 4))], [], []])
    contacts = iter([(2.0, (0, 9)), (0.5, (0, 9)), (2.0, (0, 9))])
    monkeypatch.setattr(polymer_generator_module, "_mbuild_chain", fake_chain)
    monkeypatch.setattr(polymer_generator_module, "_map_coordinates", fake_coordinates)
    monkeypatch.setattr(polymer_generator_module, "relax_in_vacuum", lambda molecule, ff: None)
    monkeypatch.setattr(polymer_generator_module, "find_pierced_rings", lambda m, x: next(checks))
    monkeypatch.setattr(polymer_generator_module, "closest_contact", lambda m, x: next(contacts))
    generator = _generator(tmp_path, max_retries=3)
    monkeypatch.setattr(generator, "_charger", Charger)

    molecule = generator.build_molecule("HABH")

    assert len(seeds) == 3 and seeds == list(range(seeds[0], seeds[0] + 3))
    assert {atom.metadata["residue_name"] for atom in molecule.atoms} == {"MMA", "HEM"}
    # The hydrogen caps have no residue_name, so they join their end units: 4 - 2 residues.
    assert {atom.metadata["residue_number"] for atom in molecule.atoms} == {1, 2}

    seeds.clear()
    checks = iter([[((0, 1, 2), (3, 4))]] * 2)
    contacts = iter([(2.0, (0, 9))] * 2)
    with pytest.raises(PolymerGenerationError, match="no build out of 2") as error:
        _generator(tmp_path, max_retries=2).build_molecule("HABH")
    assert "pierces ring" in str(error.value) and "fix:" in str(error.value)
    assert MIN_CONTACT_ANGSTROM == 1.0


def test_a_failed_build_is_reported_as_a_generation_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    generator = _generator(tmp_path)

    def fail(sequence):
        raise RuntimeError("boom")

    monkeypatch.setattr(generator, "build_molecule", fail)
    with pytest.raises(
        PolymerGenerationError, match="'HABH' could not be built: RuntimeError: boom"
    ):
        generator.generate_polymer("HABH")
    assert not generator.sdf_path("HABH").exists()


def _stub_build(monkeypatch: pytest.MonkeyPatch, generator: PolymerGenerator) -> None:
    """Replace mBuild, the charger and the relaxation so build_molecule runs in milliseconds."""

    class Charger:
        def charge_molecule(self, molecule):
            return molecule

    monkeypatch.setattr(polymer_generator_module, "_mbuild_chain", lambda *args: (None, None, None))
    monkeypatch.setattr(
        polymer_generator_module,
        "_map_coordinates",
        lambda chain, *args: np.random.default_rng(0).random((chain.GetNumAtoms(), 3)) * 10,
    )
    monkeypatch.setattr(polymer_generator_module, "relax_in_vacuum", lambda molecule, ff: None)
    monkeypatch.setattr(polymer_generator_module, "find_pierced_rings", lambda m, x: [])
    monkeypatch.setattr(polymer_generator_module, "closest_contact", lambda m, x: (2.0, (0, 1)))
    monkeypatch.setattr(generator, "_charger", Charger)


def _residue_names(molecule) -> list[tuple[int, str, str]]:
    """Return the (number, residue_name, extended_name) of each residue in atom order."""
    seen = []
    for atom in molecule.atoms:
        key = (
            atom.metadata["residue_number"],
            atom.metadata["residue_name"],
            atom.metadata["extended_name"],
        )
        if not seen or seen[-1] != key:
            seen.append(key)
    return seen


def test_caps_with_a_residue_name_are_residues_of_their_own(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A 5-mer with two own-residue caps has 5 residues; the caps are first and last."""
    generator = PolymerGenerator.from_smiles(
        [("A", "MMA", "MMA", MMA), ("B", "HEMA", "HEM", HEMA)],
        head=[_option("M", METHYL_CAP, residue_name="MET", name="methyl")],
        tail=[_option("V", "[*]C=C(C)C(=O)OC", residue_name="MMV")],
        cache_directory=tmp_path,
        force_field=SAGE,
    )
    _stub_build(monkeypatch, generator)

    molecule = generator.build_molecule("MABAV")

    assert _residue_names(molecule) == [
        (1, "MET", "methyl"),
        (2, "MMA", "MMA"),
        (3, "HEM", "HEMA"),
        (4, "MMA", "MMA"),
        (5, "MMV", "V"),  # a cap without a name is named by its label
    ]


def test_caps_without_a_residue_name_join_their_end_units(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A 5-mer with two folded caps has 3 residues, the repeat units."""
    generator = PolymerGenerator.from_smiles(
        [("A", "MMA", "MMA", MMA), ("B", "HEMA", "HEM", HEMA)],
        head=[_option("M", METHYL_CAP)],
        tail=H_OPTIONS,
        cache_directory=tmp_path,
        force_field=SAGE,
    )
    _stub_build(monkeypatch, generator)

    molecule = generator.build_molecule("MABAH")

    assert [(n, name) for n, name, _ in _residue_names(molecule)] == [
        (1, "MMA"),
        (2, "HEM"),
        (3, "MMA"),
    ]
    sizes = [
        sum(atom.metadata["residue_number"] == n for atom in molecule.atoms) for n in (1, 2, 3)
    ]
    assert sizes == [15 + 4, 19, 15 + 1]


# ---------------------------------------------------------------------------
# Sequences and the builder
# ---------------------------------------------------------------------------


def _builder(tmp_path: Path, length: int = 5, end_groups: dict | None = None) -> PolymerBuilder:
    return PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.5, 0.5],
        length=length,
        type_prefix="MH",
        cache_directory=tmp_path,
        generation_mode="dynamic",
        monomer_smiles={"MMA": MMA, "HEMA": HEMA},
        monomer_names={"A": "MMA", "B": "HEMA"},
        residue_names={"MMA": "MMA", "HEMA": "HEM"},
        end_groups=end_groups or {"head": H_OPTIONS, "tail": H_OPTIONS},
        force_field=SAGE,
    )


def _count_generations(monkeypatch: pytest.MonkeyPatch, builder: PolymerBuilder) -> list[str]:
    """Replace the generator's build with a stub that records each sequence."""
    builder._ensure_generators_initialized()
    calls: list[str] = []

    def fake_generate(sequence: str) -> str:
        calls.append(sequence)
        return f"molecule {sequence}"

    monkeypatch.setattr(builder._polymer_generator, "generate_polymer", fake_generate)
    return calls


def test_the_same_seed_draws_the_same_chains(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    first, second = _builder(tmp_path), _builder(tmp_path)
    _count_generations(monkeypatch, first)
    _count_generations(monkeypatch, second)

    first.build(count=20, seed=7)
    second.build(count=20, seed=7)

    assert first._generated_sequences == second._generated_sequences
    assert first.sequence_counts == second.sequence_counts


def test_dynamic_mode_builds_each_drawn_sequence_once_in_its_direction(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A star-sited chain has a direction, so dynamic mode keeps the drawn sequence."""
    builder = _builder(tmp_path)
    calls = _count_generations(monkeypatch, builder)

    random_state = __import__("random").Random(3)

    def draw() -> str:
        """Head cap, the 3 repeat units of a 5-mer, tail cap, in the builder's order."""
        head = random_state.choices(["H"], weights=[1.0], k=1)[0]
        middle = "".join(random_state.choices(["A", "B"], weights=[0.5, 0.5], k=3))
        return head + middle + random_state.choices(["H"], weights=[1.0], k=1)[0]

    drawn = [draw() for _ in range(30)]
    molecules, counts = builder.build(count=30, seed=3)

    sequences = builder._generated_sequences
    assert sequences == drawn
    assert any(sequence != canonical_sequence(sequence) for sequence in sequences)
    assert sorted(calls) == sorted(set(sequences))
    assert molecules == [f"molecule {sequence}" for sequence in builder.sequence_counts]
    assert sum(counts) == 30


def test_each_chain_draws_its_caps_from_their_probabilities(tmp_path: Path) -> None:
    """The head cap follows its probabilities; a fixed tail cap is on every chain."""
    import random

    end_groups = {
        "head": [_option("H", H_CAP, 0.3), _option("M", METHYL_CAP, 0.7)],
        "tail": [_option("T", METHYL_CAP, 1.0)],
    }
    builder = _builder(tmp_path, length=6, end_groups=end_groups)
    random.seed(2024)
    sequences = [builder._draw_sequence() for _ in range(2000)]

    heads = [sequence[0] for sequence in sequences]
    assert set(heads) == {"H", "M"}
    assert heads.count("H") / 2000 == pytest.approx(0.3, abs=0.03)
    assert all(sequence[-1] == "T" for sequence in sequences)
    # The length counts the two ends: a 6-mer has 4 repeat units.
    assert all(len(sequence) == 6 for sequence in sequences)
    assert all(set(sequence[1:-1]) <= {"A", "B"} for sequence in sequences)
    middles = "".join(sequence[1:-1] for sequence in sequences)
    assert middles.count("A") / len(middles) == pytest.approx(0.5, abs=0.03)


def test_built_chains_are_head_plus_repeat_units_plus_tail(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The sequence string of each built chain is its head label, units and tail label."""
    end_groups = {
        "head": [_option("H", H_CAP, 0.5), _option("M", METHYL_CAP, 0.5)],
        "tail": [_option("T", METHYL_CAP, 1.0)],
    }
    builder = _builder(tmp_path, length=5, end_groups=end_groups)
    calls = _count_generations(monkeypatch, builder)

    builder.build(count=50, seed=1)

    for sequence in builder._generated_sequences:
        head, middle, tail = builder._polymer_generator.split(sequence)
        assert head.label + middle + tail.label == sequence
        assert len(middle) == 5 - 2
        assert tail.label == "T"
    assert {sequence[0] for sequence in builder._generated_sequences} == {"H", "M"}
    assert set(calls) == set(builder._generated_sequences)
    assert builder._get_sdf_path(calls[0], tmp_path).name == f"MH_seq={calls[0]}_5-mer_charged.sdf"


def test_a_sequence_and_its_reverse_are_different_chains() -> None:
    """The head cap sits on the first unit's first '*', so AAB is not BAA."""
    units = _units(A=MMA, B=HEMA)
    forward, _ = assemble_chain("AAB", units, _cap(H_CAP), _cap(H_CAP))
    backward, _ = assemble_chain("BAA", units, _cap(H_CAP), _cap(H_CAP))
    assert Chem.MolToSmiles(forward) != Chem.MolToSmiles(backward)


def test_cached_mode_reads_a_dynamic_chain_named_in_its_drawn_direction(tmp_path: Path) -> None:
    """Cached mode looks for the canonical name, then for the reversed one."""
    reversed_file = tmp_path / "MH_seq=BAA_3-mer_charged.sdf"
    reversed_file.write_text("chain", encoding="utf-8")
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.5, 0.5],
        length=3,
        type_prefix="MH",
        sdf_directory=tmp_path,
    )
    assert builder._get_sdf_path("AAB", tmp_path) == reversed_file
    canonical_file = tmp_path / "MH_seq=AAB_3-mer_charged.sdf"
    canonical_file.write_text("chain", encoding="utf-8")
    assert builder._get_sdf_path("AAB", tmp_path) == canonical_file


def test_canonical_sequence_reads_a_chain_from_its_smaller_end() -> None:
    assert canonical_sequence("ABAAA") == "AAABA"
    assert canonical_sequence("AABBA") == "AABBA"
    assert canonical_sequence("BA") == "AB"


# ---------------------------------------------------------------------------
# Provenance and citation
# ---------------------------------------------------------------------------


def test_methods_to_cite_names_mbuild_only_when_it_built_chains() -> None:
    from types import SimpleNamespace

    from polyzymd.simulation.artifact_integrity import methods_to_cite

    config = SimpleNamespace(solvent=None)
    assert MBUILD not in methods_to_cite(config, {"polymer_sequence_seed": 1})
    assert MBUILD not in methods_to_cite(config)
    cited = methods_to_cite(config, {"polymer_builder": "mbuild 1.2.0"})
    assert MBUILD in cited
    assert "doi:10.1007/978-981-10-1128-3_5" in MBUILD


def test_mbuild_provenance_names_the_installed_release() -> None:
    import mbuild

    record = mbuild_provenance()
    assert record["package"] == "mbuild"
    assert record["version"] == mbuild.__version__
    if record["commit"] is None:
        pytest.skip("mBuild is not installed from git")
    assert record["commit"] == MBUILD_COMMIT
    assert record["requested_revision"] == "biopolymers-v0.1.0"
    assert MBUILD_COMMIT in polymer_generator_module.mbuild_summary()


# ---------------------------------------------------------------------------
# Real builds
# ---------------------------------------------------------------------------


def _checked_build(generator: PolymerGenerator, sequence: str, monkeypatch) -> object:
    """Build ``sequence``, check the files and the geometry, and read it back from the cache."""
    from openff.units import unit

    molecule = generator.generate_polymer(sequence)

    sdf = generator.sdf_path(sequence)
    assert sdf.exists()
    assert json.loads(generator.metadata_path(sdf).read_text()) == generator.cache_metadata(
        sequence
    )
    charges = molecule.partial_charges.m_as(unit.elementary_charge)
    assert abs(charges.sum() - molecule.total_charge.m_as(unit.elementary_charge)) < 1e-6
    assert abs(charges.sum()) < 1e-6

    rdmol = molecule.to_rdkit()
    xyz = molecule.conformers[0].m_as(unit.angstrom)
    assert find_pierced_rings(rdmol, xyz) == []
    assert closest_contact(rdmol, xyz)[0] >= MIN_CONTACT_ANGSTROM

    def fail(sequence):
        pytest.fail("the cached chain was built again")

    monkeypatch.setattr(generator, "build_molecule", fail)
    cached = generator.generate_polymer(sequence)
    assert np.allclose(cached.partial_charges.m_as(unit.elementary_charge), charges)
    return molecule


def _residues(molecule) -> list[tuple[int, str]]:
    """Return the (number, name) of each residue in atom order."""
    seen = []
    for atom in molecule.atoms:
        key = (atom.metadata["residue_number"], atom.metadata["residue_name"])
        if not seen or seen[-1] != key:
            seen.append(key)
    return seen


def test_a_real_mma_hema_chain_is_built_charged_and_cached(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """ABABA of MMA and HEMA builds with mBuild, NAGL and Sage, and is reused from its cache."""
    pytest.importorskip("mbuild")
    generator = _generator(tmp_path)

    molecule = _checked_build(generator, "HABABAH", monkeypatch)

    assert molecule.n_atoms == 85
    assert _residues(molecule) == [(1, "MMA"), (2, "HEM"), (3, "MMA"), (4, "HEM"), (5, "MMA")]


def test_a_real_zwitterionic_chain_keeps_its_formal_charges(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """AABAB of SBMA and EGPMA keeps three N+ and three O- through the SDF."""
    pytest.importorskip("mbuild")
    generator = PolymerGenerator.from_smiles(
        [("A", "SBMA", "SBM", SBMA), ("B", "EGPMA", "EGP", EGPMA)],
        head=[_option("M", METHYL_CAP)],
        tail=H_OPTIONS,
        cache_directory=tmp_path,
        force_field=SAGE,
    )

    molecule = _checked_build(generator, "MAABABH", monkeypatch)

    assert molecule.n_atoms == 180
    formal = [
        (atom.symbol, atom.formal_charge.m) for atom in molecule.atoms if atom.formal_charge.m
    ]
    assert sorted(formal) == [("N", 1)] * 3 + [("O", -1)] * 3
    assert [name for _, name in _residues(molecule)] == ["SBM", "SBM", "EGP", "SBM", "EGP"]


def test_a_real_methoxyethyl_methacrylate_chain_builds(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    pytest.importorskip("mbuild")
    generator = PolymerGenerator.from_smiles(
        [("A", "MEMA", "MEM", MEMA)],
        head=H_OPTIONS,
        tail=H_OPTIONS,
        cache_directory=tmp_path,
        force_field=SAGE,
    )

    molecule = _checked_build(generator, "HAAAAAH", monkeypatch)

    assert molecule.n_atoms == 5 * 22 + 2


def test_a_real_sbma_egma_chain_with_vinyl_caps_builds(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The Paper 1 chemistry: vinyl monomer caps, each its own residue named after its monomer.

    AABAB is head cap A (SBM vinyl), repeat units A, B, A and tail cap B (EGM vinyl).
    """
    pytest.importorskip("mbuild")
    generator = PolymerGenerator.from_smiles(
        [("A", "SBMA", "SBM", SBM), ("B", "EGMA", "EGM", EGM)],
        head=[
            _option("A", SBM_VINYL, 0.5, residue_name="SBM"),
            _option("B", EGM_VINYL, 0.5, residue_name="EGM"),
        ],
        tail=[
            _option("A", SBM_VINYL, 0.5, residue_name="SBM"),
            _option("B", EGM_VINYL, 0.5, residue_name="EGM"),
        ],
        cache_directory=tmp_path,
        force_field=SAGE,
        type_prefix="SBMA-EGMA",
    )

    molecule = _checked_build(generator, "AABAB", monkeypatch)

    assert generator.sdf_path("AABAB").name == "SBMA-EGMA_seq=AABAB_5-mer_charged.sdf"
    # Caps: SBM vinyl 38, EGM vinyl 42; units: SBM 39, EGM 43, SBM 39.
    assert molecule.n_atoms == 38 + 39 + 43 + 39 + 42 == 201
    assert [name for _, name in _residues(molecule)] == ["SBM", "SBM", "EGM", "SBM", "EGM"]
    assert [number for number, _ in _residues(molecule)] == [1, 2, 3, 4, 5]
    assert sum(atom.formal_charge.m for atom in molecule.atoms) == 0
    formal = sorted(
        (atom.symbol, atom.formal_charge.m) for atom in molecule.atoms if atom.formal_charge.m
    )
    assert formal == [("N", 1)] * 3 + [("O", -1)] * 3
    rdmol = Chem.RemoveHs(molecule.to_rdkit())
    assert len(rdmol.GetSubstructMatches(Chem.MolFromSmarts("[CH1]=[CX3](C)C(=O)O"))) == 2


def _backbone_handedness(sequence: str, monomers: dict[str, MonomerUnit]) -> list[int]:
    """Return the sign of each repeat unit's backbone quaternary carbon in an mBuild chain.

    The sign is that of the triple product of the bonds from the quaternary
    carbon to its own CH2, to the next unit's CH2 (the tail cap on the last
    unit) and to its ester carbon. Chlorine caps keep each unit's CH2 a CH2.
    """
    head, tail = _cap("[*]Cl"), _cap("[*]Cl")
    chain, residues = assemble_chain(sequence, monomers, head, tail)
    elements, bonds, xyz = polymer_generator_module._mbuild_chain(sequence, monomers, head, tail, 7)
    xyz = polymer_generator_module._map_coordinates(chain, elements, bonds, xyz)

    def hydrogens(atom) -> int:
        return atom.GetTotalNumHs(includeNeighbors=True)

    signs = {}
    for atom in chain.GetAtoms():
        if atom.GetSymbol() != "C" or atom.GetDegree() != 4 or hydrogens(atom):
            continue
        quaternary, residue = atom.GetIdx(), residues[atom.GetIdx()]
        neighbors = list(atom.GetNeighbors())
        ester = [
            n.GetIdx()
            for n in neighbors
            if any(
                b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol() == "O"
                for b in n.GetBonds()
            )
        ]
        own = [
            n.GetIdx()
            for n in neighbors
            if n.GetSymbol() == "C" and hydrogens(n) == 2 and residues[n.GetIdx()] == residue
        ]
        methyl = [n.GetIdx() for n in neighbors if n.GetSymbol() == "C" and hydrogens(n) == 3]
        assert len(ester) == len(own) == len(methyl) == 1
        [following] = [n.GetIdx() for n in neighbors if n.GetIdx() not in ester + own + methyl]
        vectors = [xyz[i] - xyz[quaternary] for i in (own[0], following, ester[0])]
        signs[residue] = int(np.sign(np.dot(vectors[0], np.cross(vectors[1], vectors[2]))))
    return [signs[residue] for residue in sorted(signs)]


def test_stereo_letters_set_the_backbone_handedness() -> None:
    """One stereo letter gives one handedness on every unit; alternating letters alternate it."""
    pytest.importorskip("mbuild")
    monomers = {
        "A": MonomerUnit("A", "MMA", "MMA", parse_star_smiles(MMA_R, 2, "monomer 'A'")),
        "B": MonomerUnit("B", "MMA", "MMA", parse_star_smiles(MMA_S, 2, "monomer 'B'")),
    }

    isotactic = _backbone_handedness("AAAAA", monomers)
    syndiotactic = _backbone_handedness("ABABA", monomers)

    assert len(isotactic) == len(syndiotactic) == 5
    assert len(set(isotactic)) == 1 and 0 not in isotactic
    assert syndiotactic == [syndiotactic[0], -syndiotactic[0]] * 2 + [syndiotactic[0]]
    assert syndiotactic[0] == isotactic[0]


@pytest.mark.parametrize(
    ("cap", "expected"), [("[*]/C=C(/C)C(=O)OCCO", "Z"), ("[*]/C=C(\\C)C(=O)OCCO", "E")]
)
def test_a_double_bond_written_in_a_cap_keeps_its_configuration(
    tmp_path: Path, cap, expected
) -> None:
    """E/Z written next to a '*' survives the build and is recorded in the SDF."""
    pytest.importorskip("mbuild")
    from rdkit.Chem import rdCIPLabeler

    option = {"label": "A", "name": "HEMA", "residue_name": "HMA", "smiles": cap}
    generator = PolymerGenerator.from_smiles(
        [("A", "HEMA", "HMA", HEMA)],
        head=[option],
        tail=[option],
        cache_directory=tmp_path,
        force_field=SAGE,
        type_prefix="T",
    )
    rdmol = generator.generate_polymer("AAAAA").to_rdkit()
    rdCIPLabeler.AssignCIPLabels(rdmol)
    codes = [
        bond.GetPropsAsDict().get("_CIPCode")
        for bond in rdmol.GetBonds()
        if bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetBeginAtom().GetSymbol() == "C"
        and bond.GetEndAtom().GetSymbol() == "C"
    ]
    assert codes == [expected, expected]
