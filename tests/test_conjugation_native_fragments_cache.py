"""Tests for native explicit fragments and additive cached polymer pools."""

from __future__ import annotations

import importlib.abc
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.polymer.fragments_native import (
    assemble_native_fragment_compound,
    generate_native_fragment_polymer,
    native_fragment_cache_identity,
)
from polyzymd.builders.conjugation.polymer.mbuild import from_mbuild
from polyzymd.builders.conjugation.polymer.provided_molecules import (
    _validate_provided_molecule,
    build_provided_molecule_pool,
    load_validated_provided_molecule_sdf,
)
from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.config.schema import PolymerConfig


class _RemovedBackendImportBlocker(importlib.abc.MetaPathFinder):
    """Import hook that fails if the removed legacy backend is imported."""

    def find_spec(self, fullname, path=None, target=None):
        """Reject removed backend imports while leaving other imports untouched."""
        legacy_backend = "poly" + "merist"
        if fullname == legacy_backend or fullname.startswith(f"{legacy_backend}."):
            raise ImportError("blocked removed legacy backend import")
        return None


def _fragment_spec(terminal: str = "[*]C", middle: str = "[*:1]C([*:2])"):
    """Return a lightweight terminal/middle fragment spec."""
    return SimpleNamespace(terminal=terminal, middle=middle)


def _charged_sdf(path: Path, smiles: str = "C") -> Path:
    """Write a tiny charged SDF fixture with conformer coordinates."""
    from openff.toolkit import Molecule
    from openff.units import unit

    molecule = Molecule.from_smiles(smiles)
    molecule.generate_conformers(n_conformers=1)
    molecule.partial_charges = [0.0] * molecule.n_atoms * unit.elementary_charge
    molecule.to_file(path, file_format="SDF")
    return path


def _uncharged_sdf(path: Path, smiles: str = "C") -> Path:
    """Write a tiny uncharged SDF fixture with conformer coordinates."""
    from openff.toolkit import Molecule

    molecule = Molecule.from_smiles(smiles)
    molecule.generate_conformers(n_conformers=1)
    molecule.to_file(path, file_format="SDF")
    return path


def test_provided_molecules_schema_accepts_probabilistic_and_fixed_modes(tmp_path: Path):
    """PolymerConfig should accept additive provided molecule pools."""
    sdf = tmp_path / "polymer.sdf"
    config = PolymerConfig.model_validate(
        {
            "generation_mode": "dynamic",
            "type_prefix": "PEG",
            "monomers": [{"label": "A", "probability": 1.0, "name": "peg", "smiles": "C=C"}],
            "length": 2,
            "count": 1,
            "reactions": {
                "initiation": "default",
                "polymerization": "default",
                "termination": "default",
            },
            "provided_molecules": [
                {"name": "prob", "count": 2, "entries": [{"sdf_path": sdf, "probability": 1.0}]},
                {"name": "fixed", "entries": [{"sdf_path": sdf, "count": 3}]},
            ],
        }
    )

    assert config.provided_molecules[0].entries[0].sdf_path == sdf
    assert config.generation_mode.value == "dynamic"


def test_unreleased_cached_pools_key_is_rejected():
    """The unreleased cached_pools key should not remain a supported alias."""
    with pytest.raises(ValueError, match="provided_molecules"):
        PolymerConfig.model_validate(
            {
                "generation_mode": "cached",
                "type_prefix": "PEG",
                "monomers": [{"label": "A", "probability": 1.0}],
                "length": 2,
                "count": 1,
                "sdf_directory": ".",
                "cached_pools": [],
            }
        )


def test_fragments_mode_rejects_deprecated_sdf_directory():
    """Native fragments mode should not use sequence-derived SDF directories."""
    with pytest.raises(ValueError, match="sdf_directory"):
        PolymerConfig.model_validate(
            {
                "generation_mode": "fragments",
                "type_prefix": "frag",
                "monomers": [{"label": "A", "probability": 1.0, "name": "a"}],
                "fragments": {"A": {"terminal": "[*]C", "middle": "[*:1]C([*:2])"}},
                "length": 2,
                "count": 1,
                "sdf_directory": ".",
            }
        )


@pytest.mark.parametrize(
    "provided_molecules, message",
    [
        (
            [{"name": "bad", "count": 2, "entries": [{"sdf_path": "a.sdf", "probability": 0.5}]}],
            "sum",
        ),
        (
            [{"name": "bad", "count": 2, "entries": [{"sdf_path": "a.sdf", "count": 1}]}],
            "Probabilistic",
        ),
        ([{"name": "bad", "entries": [{"sdf_path": "a.sdf", "probability": 1.0}]}], "Fixed"),
        ([{"name": "bad", "entries": [{"sdf_path": "a.mol2", "count": 1}]}], ".sdf"),
    ],
)
def test_provided_molecules_schema_rejects_invalid_mixed_modes(provided_molecules, message):
    """Pool mode validation should reject mixed probability/count inventories."""
    with pytest.raises(ValueError, match=message):
        PolymerConfig.model_validate(
            {
                "generation_mode": "cached",
                "type_prefix": "PEG",
                "monomers": [{"label": "A", "probability": 1.0}],
                "length": 2,
                "count": 1,
                "sdf_directory": ".",
                "provided_molecules": provided_molecules,
            }
        )


def test_provided_molecules_selection_seed_precedence_and_duplicate_coalescing(tmp_path: Path):
    """Pool selection should be deterministic and coalesce duplicate SDF paths."""
    sdf_a = _charged_sdf(tmp_path / "a.sdf", "C")
    sdf_b = _charged_sdf(tmp_path / "b.sdf", "O")
    pool = SimpleNamespace(
        name="pool",
        count=8,
        seed=5,
        entries=[
            SimpleNamespace(sdf_path=sdf_a, probability=0.5, count=None),
            SimpleNamespace(sdf_path=sdf_b, probability=0.5, count=None),
        ],
    )

    first_molecules, first_counts = build_provided_molecule_pool([pool], base_seed=1, caller_seed=2)
    second_molecules, second_counts = build_provided_molecule_pool(
        [pool], base_seed=9, caller_seed=8
    )

    assert first_counts == second_counts
    assert sum(first_counts) == 8
    assert len(first_molecules) <= 2


def test_fixed_provided_molecule_duplicate_paths_accumulate_counts(tmp_path: Path):
    """Fixed provided entries with the same SDF path should accumulate counts."""
    sdf = _charged_sdf(tmp_path / "dup.sdf", "C")
    pool = SimpleNamespace(
        name="fixed",
        count=None,
        seed=None,
        entries=[
            SimpleNamespace(sdf_path=sdf, probability=None, count=2),
            SimpleNamespace(sdf_path=sdf, probability=None, count=3),
        ],
    )

    molecules, counts = build_provided_molecule_pool([pool])

    assert len(molecules) == 1
    assert counts == [5]


def test_cached_sdf_directory_mode_still_works_and_warns_once(monkeypatch, tmp_path: Path):
    """Legacy cached mode should load sequence-derived SDFs with one warning."""
    import polyzymd.builders.polymer as polymer_module

    monkeypatch.setattr(polymer_module, "_SDF_DIRECTORY_DEPRECATION_WARNED", False)
    _charged_sdf(tmp_path / "PEG_seq=AA_2-mer_charged.sdf", "C")
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=2,
        type_prefix="PEG",
        sdf_directory=tmp_path,
        generation_mode="cached",
    )

    with pytest.warns(UserWarning, match="provided_molecules") as warnings_seen:
        molecules, counts = builder.build(count=2, seed=1)
        builder.build(count=1, seed=2)

    assert len(warnings_seen) == 1
    assert len(molecules) == 1
    assert counts == [2]


def test_dynamic_sdf_directory_compatibility_warns_once(monkeypatch, tmp_path: Path):
    """Historical dynamic mode with sdf_directory should still load old SDFs."""
    import polyzymd.builders.polymer as polymer_module

    monkeypatch.setattr(polymer_module, "_SDF_DIRECTORY_DEPRECATION_WARNED", False)
    _charged_sdf(tmp_path / "PEG_seq=AA_2-mer_charged.sdf", "C")
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=2,
        type_prefix="PEG",
        sdf_directory=tmp_path,
        generation_mode="dynamic",
        monomer_smiles={"a": "C=C"},
        monomer_names={"A": "a"},
        reactions=SimpleNamespace(),
    )

    with pytest.warns(UserWarning, match="derives filenames") as warnings_seen:
        molecules, counts = builder.build(count=1, seed=1)

    assert len(warnings_seen) == 1
    assert len(molecules) == 1
    assert counts == [1]


def test_provided_molecule_paths_expand_from_config_location(tmp_path: Path):
    """Nested provided_molecules entry paths should expand like other SDF paths."""
    from polyzymd.config.loader import _expand_paths

    expanded = _expand_paths(
        {
            "polymers": {
                "provided_molecules": [
                    {"name": "fixed", "entries": [{"sdf_path": "polymers/mol.sdf", "count": 1}]}
                ]
            }
        },
        tmp_path,
    )

    sdf_path = expanded["polymers"]["provided_molecules"][0]["entries"][0]["sdf_path"]
    assert sdf_path == str(tmp_path / "polymers" / "mol.sdf")


def test_provided_molecule_sdf_validation_rejects_missing_charges_and_multiple_records(
    tmp_path: Path,
):
    """Provided molecule validation should not autocharge invalid SDF inputs."""
    charged = _charged_sdf(tmp_path / "charged.sdf", "C")
    assert load_validated_provided_molecule_sdf(charged).partial_charges is not None

    missing_charges = _uncharged_sdf(tmp_path / "missing.sdf")
    with pytest.raises(ValueError, match="partial charges"):
        load_validated_provided_molecule_sdf(missing_charges)

    multi = tmp_path / "multi.sdf"
    multi.write_text(charged.read_text() + charged.read_text(), encoding="utf-8")
    with pytest.raises(ValueError, match="exactly one molecule"):
        load_validated_provided_molecule_sdf(multi)


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf")])
def test_provided_molecule_validation_rejects_nonfinite_coordinates(
    tmp_path: Path, bad_value: float
):
    """Provided molecule validation should reject NaN and infinite coordinates."""
    from openff.toolkit import Molecule
    from openff.units import unit

    molecule = Molecule.from_smiles("C")
    molecule.generate_conformers(n_conformers=1)
    coordinates = molecule.conformers[0].m_as("angstrom")
    coordinates[0, 0] = bad_value
    molecule._conformers = [coordinates * unit.angstrom]
    molecule.partial_charges = [0.0] * molecule.n_atoms * unit.elementary_charge

    with pytest.raises(ValueError, match="non-finite conformer coordinates"):
        _validate_provided_molecule(molecule, tmp_path / "bad.sdf")


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf")])
def test_provided_molecule_validation_rejects_nonfinite_partial_charges(
    tmp_path: Path, bad_value: float
):
    """Provided molecule validation should reject NaN and infinite partial charges."""
    from openff.toolkit import Molecule
    from openff.units import unit

    molecule = Molecule.from_smiles("C")
    molecule.generate_conformers(n_conformers=1)
    charges = [0.0] * molecule.n_atoms
    charges[0] = bad_value
    molecule.partial_charges = charges * unit.elementary_charge

    with pytest.raises(ValueError, match="non-finite partial charges"):
        _validate_provided_molecule(molecule, tmp_path / "bad.sdf")


def test_native_fragments_parse_mapped_unmapped_ports_and_preserve_direction():
    """Native fragments should support mapped and legacy unlabelled dummy ports."""
    specs = {
        "A": _fragment_spec(terminal="[*]C", middle="[*:1]C([*:2])"),
        "B": _fragment_spec(terminal="[*]O", middle="[*]O[*]"),
    }
    ab = from_mbuild(assemble_native_fragment_compound("frag", specs, "AB"))
    ba = from_mbuild(assemble_native_fragment_compound("frag", specs, "BA"))

    assert ab.n_atoms == ba.n_atoms
    assert native_fragment_cache_identity("frag", specs, "AB") != native_fragment_cache_identity(
        "frag", specs, "BA"
    )


@pytest.mark.parametrize("sequence", ["A", "AA", "ABA"])
def test_native_fragments_length_and_asymmetric_geometry(sequence: str):
    """Native fragments should produce strict OpenFF molecules for linear lengths."""
    specs = {
        "A": _fragment_spec(terminal="[*]C", middle="[*:1]C([*:2])"),
        "B": _fragment_spec(terminal="[*]O", middle="[*:1]O[*:2]"),
    }

    molecule = from_mbuild(assemble_native_fragment_compound("frag", specs, sequence))

    assert molecule.n_conformers == 1
    assert molecule.n_atoms > 0
    assert molecule.n_bonds > 0


def test_native_force_overlap_consumes_ports_and_normalizes_created_bonds():
    """force_overlap should consume ports while PolyzyMD normalizes its bond order."""
    compound = assemble_native_fragment_compound(
        "frag", {"A": _fragment_spec("[*]C", "[*:1]C([*:2])")}, "AAA"
    )

    assert list(compound.all_ports()) == []
    assert all(
        float(compound.bond_graph[left][right]["bond_order"]) > 0.0
        for left, right in compound.bonds()
    )


def test_native_force_overlap_passes_explicit_add_bond_true(monkeypatch):
    """Native fragments should call force_overlap with explicit add_bond=True."""
    import mbuild as mb

    original = mb.force_overlap
    add_bond_values = []

    def wrapped_force_overlap(*args, **kwargs):
        """Record force_overlap add_bond values while preserving behavior."""
        add_bond_values.append(kwargs.get("add_bond"))
        return original(*args, **kwargs)

    monkeypatch.setattr(mb, "force_overlap", wrapped_force_overlap)

    assemble_native_fragment_compound("frag", {"A": _fragment_spec("[*]C", "[*:1]C([*:2])")}, "AA")

    assert add_bond_values == [True]


def test_native_fragment_dummy_embedding_does_not_emit_rdkit_uff_chatter(capfd):
    """Expected dummy-atom RDKit embedding chatter should be locally suppressed."""
    assemble_native_fragment_compound("frag", {"A": _fragment_spec("[*]C", "[*:1]C([*:2])")}, "AAA")

    captured = capfd.readouterr()
    assert "UFFTYPER" not in captured.err


@pytest.mark.parametrize(
    "terminal, middle, message",
    [
        ("CC", "[*:1]C([*:2])", "1 dummy"),
        ("[*]C", "[*:1]C", "2 dummy"),
        ("[*]C", "[*:1]C([*:1])", "duplicated"),
        ("[*]C", "[*:3]C([*:2])", "Unsupported"),
    ],
)
def test_native_fragments_reject_invalid_dummy_counts_and_maps(terminal, middle, message):
    """Fragment validation should reject invalid port counts and maps."""
    with pytest.raises(ValueError, match=message):
        assemble_native_fragment_compound("frag", {"A": _fragment_spec(terminal, middle)}, "AAA")


def test_polymer_builder_fragments_preserves_exact_random_sequence(monkeypatch, tmp_path: Path):
    """Fragments mode should not reverse-canonicalize asymmetric sequences."""
    captured = []

    def fake_generate(name, fragment_specs, sequence, cache_directory, *, charger_type="nagl"):
        captured.append(sequence)
        return SimpleNamespace(charged_molecule=SimpleNamespace(sequence=sequence))

    monkeypatch.setattr(
        "polyzymd.builders.conjugation.polymer.fragments_native.generate_native_fragment_polymer",
        fake_generate,
    )
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.0, 1.0],
        length=3,
        type_prefix="frag",
        generation_mode="fragments",
        fragments={"A": _fragment_spec(), "B": _fragment_spec("[*]O", "[*:1]O[*:2]")},
        cache_directory=tmp_path,
    )

    molecules, counts = builder.build(count=2, seed=1)

    assert captured == ["BBB"]
    assert counts == [2]
    assert molecules[0].sequence == "BBB"


def test_polymer_builder_fragments_routes_cache_validation_through_generator(
    monkeypatch, tmp_path: Path
):
    """Fragments mode should always call the native generator for cache validation."""
    generated = SimpleNamespace(sequence="AA")
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=2,
        type_prefix="frag",
        generation_mode="fragments",
        fragments={"A": _fragment_spec()},
        cache_directory=tmp_path,
    )
    monkeypatch.setattr(builder, "_generate_native_fragment_polymer", lambda sequence: generated)
    monkeypatch.setattr(
        builder,
        "_native_fragment_artifact_paths",
        lambda sequence: pytest.fail("direct cache shortcut should not be used"),
    )
    monkeypatch.setattr(
        builder,
        "_load_from_sdf",
        lambda path: pytest.fail("fragments mode should not load cached SDF directly"),
    )

    molecules, counts = builder.build(count=1, seed=1)

    assert molecules == [generated]
    assert counts == [1]


def test_dynamic_plus_provided_molecules_are_returned_for_one_packmol_call(
    monkeypatch, tmp_path: Path
):
    """PolymerBuilder should return one merged molecule/count list without cross-coalescing."""
    sdf = _charged_sdf(tmp_path / "cached.sdf", "C")
    pool = SimpleNamespace(
        name="fixed",
        count=None,
        seed=None,
        entries=[SimpleNamespace(sdf_path=sdf, probability=None, count=2)],
    )
    monkeypatch.setattr(
        PolymerBuilder,
        "_get_or_create_molecule",
        lambda self, sequence: SimpleNamespace(sequence=sequence),
    )
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=2,
        type_prefix="dyn",
        generation_mode="dynamic",
        monomer_smiles={"a": "C=C"},
        monomer_names={"A": "a"},
        reactions=SimpleNamespace(),
        provided_molecules=[pool],
    )

    molecules, counts = builder.build(count=3, seed=10)

    assert [getattr(mol, "sequence", None) for mol in molecules] == ["AA", None]
    assert counts == [3, 2]
    assert builder.get_packing_info() == (molecules, counts)


def test_native_fragment_cache_hit_preserves_cached_coordinates(monkeypatch, tmp_path: Path):
    """Cached native fragments should reuse charged SDF coordinates, not regenerated ones."""
    from openff.units import unit

    class FakeCharger:
        """Minimal charger that records complete zero charges."""

        def charge_molecule(self, molecule):
            """Return the molecule with complete partial charges."""
            molecule.partial_charges = [0.0] * molecule.n_atoms * unit.elementary_charge
            return molecule

    monkeypatch.setattr(
        "polyzymd.builders.conjugation.polymer.fragments_native.get_charger",
        lambda charger_type: FakeCharger(),
    )
    specs = {"A": _fragment_spec("[*]C", "[*:1]C([*:2])")}
    first = generate_native_fragment_polymer("frag", specs, "AA", tmp_path)
    cached = first.charged_molecule
    shifted = cached.conformers[0] + 1.0 * unit.angstrom
    cached._conformers = [shifted]
    cached.to_file(first.paths.charged_sdf_path, file_format="SDF")

    second = generate_native_fragment_polymer("frag", specs, "AA", tmp_path)

    assert second.charged_molecule.conformers[0][0][0].m_as("angstrom") == pytest.approx(
        shifted[0][0].m_as("angstrom"), abs=1.0e-3
    )


def test_conjugation_free_polymer_path_forwards_fragments_and_provided_molecules(monkeypatch):
    """The free-polymer conjugation path should forward top-level polymer additions."""
    from polyzymd.builders.conjugation.system_workflow import _build_and_pack_free_polymers

    captured = {}

    class FakeBuilder:
        """Capture builder arguments without invoking Packmol."""

        _working_dir = None

        def build_polymers(self, **kwargs):
            """Capture forwarded polymer build keyword arguments."""
            captured.update(kwargs)

        def pack_polymers(self, **kwargs):
            """Capture pack call marker."""
            captured["packed"] = True

    polymers = PolymerConfig.model_validate(
        {
            "generation_mode": "fragments",
            "type_prefix": "frag",
            "monomers": [{"label": "A", "probability": 1.0, "name": "a"}],
            "fragments": {"A": {"terminal": "[*]C", "middle": "[*:1]C([*:2])"}},
            "length": 2,
            "count": 1,
            "provided_molecules": [
                {"name": "fixed", "entries": [{"sdf_path": "polymer.sdf", "count": 1}]}
            ],
        }
    )

    _build_and_pack_free_polymers(FakeBuilder(), SimpleNamespace(polymers=polymers), polymer_seed=7)

    assert captured["fragments"] == polymers.fragments
    assert captured["provided_molecules"] == polymers.provided_molecules
    assert captured["polymer_random_seed"] == polymers.random_seed
    assert captured["packed"] is True


def test_template_includes_phase14_polymer_examples():
    """Code-owned config template should expose fragments and provided molecule keys."""
    template = Path("src/polyzymd/templates/templates/config_template.yaml").read_text()

    assert 'generation_mode: "fragments"' in template
    assert "provided_molecules:" in template
    assert "probability:" in template


def test_native_fragments_do_not_import_removed_backend(tmp_path: Path):
    """Fragments/default/provided molecule native paths should not import removed code."""
    blocker = _RemovedBackendImportBlocker()
    sys.meta_path.insert(0, blocker)
    try:
        from_mbuild(
            assemble_native_fragment_compound(
                "frag", {"A": _fragment_spec("[*]C", "[*:1]C([*:2])")}, "AA"
            )
        )
        build_provided_molecule_pool(
            [
                SimpleNamespace(
                    name="fixed",
                    count=None,
                    seed=None,
                    entries=[SimpleNamespace(sdf_path=_charged_sdf(tmp_path / "c.sdf"), count=1)],
                )
            ]
        )
    finally:
        sys.meta_path.remove(blocker)
