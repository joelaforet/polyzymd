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
        head=H_CAP,
        tail=H_CAP,
        **arguments,
    )


def test_from_smiles_reads_the_monomers_and_caps(tmp_path: Path) -> None:
    generator = _generator(tmp_path)
    hema = generator.monomers["B"]
    assert (hema.label, hema.name, hema.residue_name) == ("B", "HEMA", "HEM")
    assert hema.unit.smiles == HEMA and hema.unit.anchors == (1, 2)
    assert generator.end_groups["head"].smiles == H_CAP
    assert generator.charger_type == "nagl"
    assert generator.nagl_model == "openff-gnn-am1bcc-1.0.0.pt"
    assert _generator(tmp_path, charger_type="AM1BCC").nagl_model is None


def test_from_smiles_names_a_bad_monomer(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="^monomer 'B'"):
        PolymerGenerator.from_smiles(
            [("A", "MMA", "MMA", MMA), ("B", "X", "XXX", "[*]CC")],
            head=H_CAP,
            tail=H_CAP,
            cache_directory=tmp_path,
            force_field=SAGE,
        )


def test_the_sdf_and_its_sidecar_are_named_after_the_sequence(tmp_path: Path) -> None:
    generator = _generator(tmp_path, type_prefix="MH")
    sdf = generator.sdf_path("ABABA")
    assert sdf == tmp_path / "MH_seq=ABABA_5-mer_charged.sdf"
    assert generator.metadata_path(sdf) == tmp_path / "MH_seq=ABABA_5-mer_charged.sdf.metadata.json"


def test_the_cache_metadata_records_what_the_chain_was_built_from(tmp_path: Path) -> None:
    metadata = _generator(tmp_path).cache_metadata("AAA")

    assert metadata["schema_version"] == DYNAMIC_CACHE_SCHEMA_VERSION == 5
    assert metadata["sequence"] == "AAA"
    assert metadata["monomers"] == {"A": {"name": "MMA", "residue_name": "MMA", "smiles": MMA}}
    assert metadata["end_groups"] == {"head": H_CAP, "tail": H_CAP}
    assert metadata["relaxation_force_field"] == SAGE
    assert metadata["builder"]["package"] == "mbuild"
    assert json.loads(json.dumps(metadata)) == metadata


def test_a_sequence_with_unknown_labels_is_refused(tmp_path: Path) -> None:
    with pytest.raises(PolymerGenerationError, match="labels \\['C'\\]"):
        _generator(tmp_path).build_molecule("ACA")


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

    molecule = generator.build_molecule("AB")

    assert len(seeds) == 3 and seeds == list(range(seeds[0], seeds[0] + 3))
    assert {atom.metadata["residue_name"] for atom in molecule.atoms} == {"MMA", "HEM"}

    seeds.clear()
    checks = iter([[((0, 1, 2), (3, 4))]] * 2)
    contacts = iter([(2.0, (0, 9))] * 2)
    with pytest.raises(PolymerGenerationError, match="no build out of 2") as error:
        _generator(tmp_path, max_retries=2).build_molecule("AB")
    assert "pierces ring" in str(error.value) and "fix:" in str(error.value)
    assert MIN_CONTACT_ANGSTROM == 1.0


def test_a_failed_build_is_reported_as_a_generation_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    generator = _generator(tmp_path)

    def fail(sequence):
        raise RuntimeError("boom")

    monkeypatch.setattr(generator, "build_molecule", fail)
    with pytest.raises(PolymerGenerationError, match="'AB' could not be built: RuntimeError: boom"):
        generator.generate_polymer("AB")
    assert not generator.sdf_path("AB").exists()


# ---------------------------------------------------------------------------
# Sequences and the builder
# ---------------------------------------------------------------------------


def _builder(tmp_path: Path, length: int = 5) -> PolymerBuilder:
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
        end_groups={"head": H_CAP, "tail": H_CAP},
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


def test_each_unique_canonical_sequence_is_generated_once(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    builder = _builder(tmp_path)
    calls = _count_generations(monkeypatch, builder)

    molecules, counts = builder.build(count=30, seed=3)

    sequences = builder._generated_sequences
    assert all(sequence == canonical_sequence(sequence) for sequence in sequences)
    assert sorted(calls) == sorted(set(sequences))
    assert molecules == [f"molecule {sequence}" for sequence in builder.sequence_counts]
    assert sum(counts) == 30


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

    molecule = _checked_build(generator, "ABABA", monkeypatch)

    assert molecule.n_atoms == 85
    assert _residues(molecule) == [(1, "MMA"), (2, "HEM"), (3, "MMA"), (4, "HEM"), (5, "MMA")]


def test_a_real_zwitterionic_chain_keeps_its_formal_charges(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """AABAB of SBMA and EGPMA keeps three N+ and three O- through the SDF."""
    pytest.importorskip("mbuild")
    generator = PolymerGenerator.from_smiles(
        [("A", "SBMA", "SBM", SBMA), ("B", "EGPMA", "EGP", EGPMA)],
        head=METHYL_CAP,
        tail=H_CAP,
        cache_directory=tmp_path,
        force_field=SAGE,
    )

    molecule = _checked_build(generator, "AABAB", monkeypatch)

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
        head=H_CAP,
        tail=H_CAP,
        cache_directory=tmp_path,
        force_field=SAGE,
    )

    molecule = _checked_build(generator, "AAAAA", monkeypatch)

    assert molecule.n_atoms == 5 * 22 + 2
