"""Regression tests for dynamic polymer generation mapping."""

import importlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

from polyzymd.builders.polymer import PolymerBuilder

REAL_POLYMERIST_AVAILABLE = importlib.util.find_spec("polymerist") is not None


def _make_stub_package(name: str) -> ModuleType:
    """Create a package-like module stub for optional chemistry dependencies."""
    module = ModuleType(name)
    module.__path__ = []
    return module


def _install_polymerist_stubs() -> None:
    """Install minimal Polymerist stubs when the optional dependency is absent."""
    if importlib.util.find_spec("polymerist") is not None:
        return

    package_names = [
        "polymerist",
        "polymerist.genutils",
        "polymerist.genutils.fileutils",
        "polymerist.mdtools",
        "polymerist.mdtools.openfftools",
        "polymerist.polymers",
        "polymerist.rdutils",
        "polymerist.rdutils.rdcoords",
    ]
    for name in package_names:
        sys.modules.setdefault(name, _make_stub_package(name))

    pathutils = ModuleType("polymerist.genutils.fileutils.pathutils")
    pathutils.assemble_path = lambda *args, **kwargs: Path(*args)
    sys.modules.setdefault("polymerist.genutils.fileutils.pathutils", pathutils)

    partition = ModuleType("polymerist.mdtools.openfftools.partition")
    partition.partition = lambda topology: True
    sys.modules.setdefault("polymerist.mdtools.openfftools.partition", partition)

    topology = ModuleType("polymerist.mdtools.openfftools.topology")
    topology.get_largest_offmol = lambda topology: object()
    topology.topology_from_sdf = lambda path: object()
    topology.topology_to_sdf = lambda path, topology: None
    sys.modules.setdefault("polymerist.mdtools.openfftools.topology", topology)

    building = ModuleType("polymerist.polymers.building")
    building.build_linear_polymer = lambda **kwargs: object()
    building.mbmol_to_openmm_pdb = lambda *args, **kwargs: None
    building.mbmol_to_rdmol = lambda chain: object()
    sys.modules.setdefault("polymerist.polymers.building", building)

    monomers = ModuleType("polymerist.polymers.monomers")
    monomers.MonomerGroup = FakeMonomerGroup
    sys.modules.setdefault("polymerist.polymers.monomers", monomers)

    piercing = ModuleType("polymerist.rdutils.rdcoords.piercing")
    piercing.summarize_ring_piercing = lambda mol: {}
    sys.modules.setdefault("polymerist.rdutils.rdcoords.piercing", piercing)


def _install_rdkit_stubs() -> None:
    """Install minimal RDKit stubs when the optional dependency is absent."""
    if importlib.util.find_spec("rdkit") is not None:
        return

    rdkit = _make_stub_package("rdkit")
    chem = ModuleType("rdkit.Chem")
    rdkit.Chem = chem
    sys.modules.setdefault("rdkit", rdkit)
    sys.modules.setdefault("rdkit.Chem", chem)


class FakeMonomerGroup:
    """Minimal MonomerGroup stand-in for Polymerist call capture."""

    def __init__(self, monomers: dict[str, object]) -> None:
        """Initialize the fake group with named fragments."""
        self.monomers = monomers
        self.term_orient: dict[str, str] = {}


_install_polymerist_stubs()
_install_rdkit_stubs()

polymer_generator_module = importlib.import_module("polyzymd.builders.polymer_generator")
PolymerGenerator = polymer_generator_module.PolymerGenerator


def _make_generator(tmp_path: Path, fragment_names: list[str]) -> PolymerGenerator:
    """Create a PolymerGenerator with fake fragment objects."""
    monomer_group = FakeMonomerGroup({name: [f"{name}:smarts"] for name in fragment_names})
    return PolymerGenerator(monomer_group=monomer_group, cache_directory=tmp_path)


def _make_list_fragment_group(overrides: dict[str, list[str]] | None = None) -> FakeMonomerGroup:
    """Create a fake MonomerGroup with Polymerist-like list fragment entries."""
    fragments = {
        "EGMA_1-site": ["[C:1]=[C:2]", "[O:3]"],
        "EGMA_2-site": ["[C:1]-[C:2]", "[O:3]"],
    }
    fragments.update(overrides or {})
    return FakeMonomerGroup(fragments)


def _require_real_polymerist() -> None:
    """Skip when Polymerist is not installed and tests use local stubs."""
    if not REAL_POLYMERIST_AVAILABLE:
        pytest.skip("Polymerist is not available")


def _real_halogen_fragments() -> dict[str, list[str]]:
    """Create valid Polymerist SMARTS fragments with PolyzyMD monomer names."""
    _require_real_polymerist()

    from polymerist.polymers.monomers.fragments import HALOGENATED_HYDROCARBON_FRAGMENTS

    return {
        "SBMA_1-site": HALOGENATED_HYDROCARBON_FRAGMENTS["chlor_term_1"],
        "SBMA_2-site": HALOGENATED_HYDROCARBON_FRAGMENTS["chlor_mid_1"],
        "PEGMA_1-site": HALOGENATED_HYDROCARBON_FRAGMENTS["brom_term_1"],
        "PEGMA_2-site": HALOGENATED_HYDROCARBON_FRAGMENTS["brom_mid_1"],
    }


def _make_real_halogen_monomer_group(overrides: dict[str, list[str]] | None = None):
    """Create a real Polymerist MonomerGroup backed by list SMARTS data."""
    _require_real_polymerist()

    from polymerist.polymers.monomers import MonomerGroup

    fragments = _real_halogen_fragments()
    fragments.update(overrides or {})
    return MonomerGroup(monomers=fragments)


def _count_chain_halogens(chain: object) -> dict[str, int]:
    """Count diagnostic halogens in a real Polymerist-built chain."""
    _require_real_polymerist()

    from polymerist.polymers.building import mbmol_to_rdmol

    mol = mbmol_to_rdmol(chain)
    return {
        "Br": sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Br"),
        "Cl": sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl"),
    }


def _patch_structure_build(monkeypatch: pytest.MonkeyPatch, calls: list[dict[str, object]]) -> None:
    """Patch Polymerist structure helpers to avoid chemistry work."""

    def fake_build_linear_polymer(**kwargs):
        """Capture the Polymerist build arguments."""
        calls.append(kwargs)
        return object()

    def fake_write_pdb(pdb_path: Path, chain: object, resname_map: dict[str, str]) -> None:
        """Create the expected PDB output without invoking OpenMM writers."""
        pdb_path.write_text("MODEL\nEND\n", encoding="utf-8")

    monkeypatch.setattr(polymer_generator_module, "_make_monomer_group", FakeMonomerGroup)
    monkeypatch.setattr(
        polymer_generator_module, "_build_linear_polymer", fake_build_linear_polymer
    )
    monkeypatch.setattr(polymer_generator_module, "_mbmol_to_rdmol", lambda chain: object())
    monkeypatch.setattr(polymer_generator_module, "_summarize_ring_piercing", lambda mol: {})
    monkeypatch.setattr(polymer_generator_module, "_mbmol_to_openmm_pdb", fake_write_pdb)


def _write_metadata(
    generator: PolymerGenerator,
    sdf_path: Path,
    sequence: str,
    monomer_names: dict[str, str],
    overrides: dict[str, object] | None = None,
) -> None:
    """Write dynamic cache metadata with optional top-level overrides."""
    metadata = generator._build_dynamic_cache_metadata(sequence, monomer_names)
    metadata.update(overrides or {})
    generator._get_dynamic_cache_metadata_path(sdf_path).write_text(
        polymer_generator_module.json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_two_monomer_b_sequence_uses_explicit_pegma_middle_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """BBBBB should pass BBB and map B to the PEGMA 2-site fragment."""
    calls: list[dict[str, object]] = []
    _patch_structure_build(monkeypatch, calls)
    generator = _make_generator(
        tmp_path,
        ["SBMA_1-site", "SBMA_2-site", "PEGMA_1-site", "PEGMA_2-site"],
    )

    generator._build_polymer_structure("BBBBB", {"A": "SBMA", "B": "PEGMA"})

    assert calls[0]["sequence"] == "BBB"
    assert calls[0]["sequence_map"] == {"B": "PEGMA_2-site"}


def test_two_monomer_b_sequence_uses_pegma_terminal_fragments(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """BBBBB should set both terminal orientations to PEGMA 1-site fragments."""
    calls: list[dict[str, object]] = []
    _patch_structure_build(monkeypatch, calls)
    generator = _make_generator(
        tmp_path,
        ["SBMA_1-site", "SBMA_2-site", "PEGMA_1-site", "PEGMA_2-site"],
    )

    generator._build_polymer_structure("BBBBB", {"A": "SBMA", "B": "PEGMA"})

    monomer_group = calls[0]["monomers"]
    assert isinstance(monomer_group, FakeMonomerGroup)
    assert monomer_group.term_orient == {"head": "PEGMA_1-site", "tail": "PEGMA_1-site"}


def test_real_polymerist_sequence_map_keeps_b_middle_on_pegma_fragment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """BBBBB should build with PEGMA middle fragments through real Polymerist."""
    _require_real_polymerist()

    from polymerist.polymers.building import build_linear_polymer

    def fake_write_pdb(pdb_path: Path, chain: object, resname_map: dict[str, str]) -> None:
        """Create the expected PDB output after the real chain build."""
        pdb_path.write_text("MODEL\nEND\n", encoding="utf-8")

    unmapped_group = _make_real_halogen_monomer_group()
    unmapped_group.term_orient = {"head": "PEGMA_1-site", "tail": "PEGMA_1-site"}
    unmapped_chain = build_linear_polymer(
        monomers=unmapped_group,
        n_monomers=5,
        sequence="BBB",
        energy_minimize=False,
        allow_partial_sequences=True,
    )
    assert _count_chain_halogens(unmapped_chain) == {"Br": 2, "Cl": 3}

    monkeypatch.setattr(polymer_generator_module, "_mbmol_to_openmm_pdb", fake_write_pdb)
    generator = PolymerGenerator(
        monomer_group=_make_real_halogen_monomer_group(),
        cache_directory=tmp_path,
        max_retries=1,
    )

    chain, pdb_path = generator._build_polymer_structure(
        "BBBBB",
        {"A": "SBMA", "B": "PEGMA"},
    )

    assert pdb_path.exists()
    assert _count_chain_halogens(chain) == {"Br": 5, "Cl": 0}


def test_single_monomer_sequence_uses_explicit_egma_middle_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """AAAAA should pass an explicit EGMA mapping instead of Polymerist defaults."""
    calls: list[dict[str, object]] = []
    _patch_structure_build(monkeypatch, calls)
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    generator._build_polymer_structure("AAAAA", {"A": "EGMA"})

    assert calls[0]["sequence"] == "AAA"
    assert calls[0]["sequence_map"] == {"A": "EGMA_2-site"}


def test_stale_dynamic_cache_without_metadata_is_not_loaded(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A dynamic charged SDF without metadata should not be directly reused."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")

    def fail_load_cache(path: Path) -> object:
        """Fail if stale SDF loading is attempted."""
        pytest.fail(f"Stale dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


@pytest.mark.parametrize("metadata", [[], "stale", 1])
def test_dynamic_cache_rejects_non_object_json_metadata(
    metadata: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache metadata sidecars must contain a JSON object."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    generator._get_dynamic_cache_metadata_path(cache_path).write_text(
        polymer_generator_module.json.dumps(metadata),
        encoding="utf-8",
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if non-object metadata cache loading is attempted."""
        pytest.fail(f"Non-object metadata dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_dynamic_cache_rejects_charger_mismatch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache metadata must include and validate the charger type."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    _write_metadata(
        generator,
        cache_path,
        "AAAAA",
        {"A": "EGMA"},
        overrides={"charger_type": "espaloma"},
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if mismatched dynamic cache loading is attempted."""
        pytest.fail(f"Mismatched dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_dynamic_cache_rejects_fragment_fingerprint_mismatch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache metadata must validate MonomerGroup fragment chemistry."""
    current_generator = PolymerGenerator(
        monomer_group=_make_list_fragment_group(),
        cache_directory=tmp_path,
    )
    stale_generator = PolymerGenerator(
        monomer_group=_make_list_fragment_group({"EGMA_2-site": ["[C:1]-[C:2]", "[N:3]"]}),
        cache_directory=tmp_path,
    )
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    _write_metadata(stale_generator, cache_path, "AAAAA", {"A": "EGMA"})

    def fail_load_cache(path: Path) -> object:
        """Fail if fingerprint-mismatched cache loading is attempted."""
        pytest.fail(f"Fingerprint-mismatched dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert current_generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_list_fragment_contents_affect_monomer_group_fingerprint(tmp_path: Path) -> None:
    """Same fragment keys with different list contents should not share a digest."""
    first_generator = PolymerGenerator(
        monomer_group=_make_list_fragment_group(),
        cache_directory=tmp_path / "first",
    )
    second_generator = PolymerGenerator(
        monomer_group=_make_list_fragment_group({"EGMA_2-site": ["[C:1]-[C:2]", "[N:3]"]}),
        cache_directory=tmp_path / "second",
    )

    first_fingerprint = first_generator._build_monomer_group_fingerprint()
    second_fingerprint = second_generator._build_monomer_group_fingerprint()

    assert first_fingerprint["algorithm"] == "polyzymd-monomer-group-sha256-v2"
    assert first_fingerprint["digest"] != second_fingerprint["digest"]


def test_real_monomer_group_list_contents_affect_fingerprint(tmp_path: Path) -> None:
    """Real Polymerist list-backed fragment content should affect the digest."""
    fragments = _real_halogen_fragments()
    first_generator = PolymerGenerator(
        monomer_group=_make_real_halogen_monomer_group(),
        cache_directory=tmp_path / "first",
    )
    second_generator = PolymerGenerator(
        monomer_group=_make_real_halogen_monomer_group({"PEGMA_2-site": fragments["SBMA_2-site"]}),
        cache_directory=tmp_path / "second",
    )

    first_fingerprint = first_generator._build_monomer_group_fingerprint()
    second_fingerprint = second_generator._build_monomer_group_fingerprint()

    assert first_fingerprint["algorithm"] == "polyzymd-monomer-group-sha256-v2"
    assert first_fingerprint["digest"] != second_fingerprint["digest"]


def test_real_monomer_group_list_content_change_rejects_cache_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache validation should reject changed real fragment lists."""
    fragments = _real_halogen_fragments()
    current_generator = PolymerGenerator(
        monomer_group=_make_real_halogen_monomer_group(),
        cache_directory=tmp_path,
    )
    stale_generator = PolymerGenerator(
        monomer_group=_make_real_halogen_monomer_group({"PEGMA_2-site": fragments["SBMA_2-site"]}),
        cache_directory=tmp_path,
    )
    cache_path = tmp_path / "PEGMA_seq=BBBBB_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    _write_metadata(stale_generator, cache_path, "BBBBB", {"A": "SBMA", "B": "PEGMA"})

    def fail_load_cache(path: Path) -> object:
        """Fail if fingerprint-mismatched cache loading is attempted."""
        pytest.fail(f"Fingerprint-mismatched dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert current_generator.get_cached_polymer("BBBBB", {"A": "SBMA", "B": "PEGMA"}) is None


def test_dynamic_generator_rejects_short_sequences(tmp_path: Path) -> None:
    """Lower-level dynamic generation should reject sequences shorter than three."""
    generator = PolymerGenerator(
        monomer_group=_make_list_fragment_group(), cache_directory=tmp_path
    )

    with pytest.raises(ValueError, match="requires sequence length >= 3"):
        generator.generate_polymer("AA", {"A": "EGMA"})


@pytest.mark.parametrize("length", [1, 2])
def test_dynamic_builder_rejects_short_lengths(length: int, tmp_path: Path) -> None:
    """Dynamic builder configuration should reject lengths below three early."""
    with pytest.raises(ValueError, match="requires length >= 3"):
        PolymerBuilder(
            characters=["A"],
            probabilities=[1.0],
            length=length,
            type_prefix="EGMA",
            cache_directory=tmp_path / "cache",
            generation_mode="dynamic",
            monomer_smiles={"EGMA": "CO"},
            monomer_names={"A": "EGMA"},
            reactions=object(),
        )


def test_dynamic_builder_routes_cache_directory_hits_through_generator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should not load cache-directory SDFs without metadata validation."""
    cache_path = tmp_path / "PEGMA_seq=BBBBB_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    sentinel = object()
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.0, 1.0],
        length=5,
        type_prefix="PEGMA",
        cache_directory=tmp_path,
        generation_mode="dynamic",
        monomer_smiles={"SBMA": "C", "PEGMA": "CO"},
        monomer_names={"A": "SBMA", "B": "PEGMA"},
        reactions=object(),
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if PolymerBuilder bypasses dynamic cache validation."""
        pytest.fail(f"Dynamic builder loaded cache directly: {path}")

    monkeypatch.setattr(builder, "_load_from_sdf", fail_load_cache)
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("BBBBB") is sentinel


def test_dynamic_builder_routes_sdf_directory_generated_artifact_through_generator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should not directly load generated-style SDF directory hits."""
    sdf_path = tmp_path / "PEGMA_seq=BBBBB_5-mer_charged.sdf"
    sdf_path.write_text("stale generated artifact", encoding="utf-8")
    sentinel = object()
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.0, 1.0],
        length=5,
        type_prefix="PEGMA",
        sdf_directory=tmp_path,
        cache_directory=tmp_path / "cache",
        generation_mode="dynamic",
        monomer_smiles={"SBMA": "C", "PEGMA": "CO"},
        monomer_names={"A": "SBMA", "B": "PEGMA"},
        reactions=object(),
    )
    builder._fragment_generator = object()
    builder._polymer_generator = PolymerGenerator(
        monomer_group=FakeMonomerGroup(
            {
                "SBMA_1-site": ["sbma-terminal"],
                "SBMA_2-site": ["sbma-middle"],
                "PEGMA_1-site": ["pegma-terminal"],
                "PEGMA_2-site": ["pegma-middle"],
            }
        ),
        cache_directory=tmp_path / "cache",
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if generated-style SDF loading bypasses validation."""
        pytest.fail(f"Dynamic builder loaded generated-style SDF directly: {path}")

    monkeypatch.setattr(builder, "_load_from_sdf", fail_load_cache)
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("BBBBB") is sentinel


def test_dynamic_builder_loads_generated_sdf_directory_hit_when_metadata_matches(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should load generated-style SDF hits with matching metadata."""
    sdf_path = tmp_path / "PEGMA_seq=BBBBB_5-mer_charged.sdf"
    sdf_path.write_text("validated generated artifact", encoding="utf-8")
    sentinel = object()
    monomer_group = FakeMonomerGroup(
        {
            "SBMA_1-site": ["sbma-terminal"],
            "SBMA_2-site": ["sbma-middle"],
            "PEGMA_1-site": ["pegma-terminal"],
            "PEGMA_2-site": ["pegma-middle"],
        }
    )
    generator = PolymerGenerator(monomer_group=monomer_group, cache_directory=tmp_path / "cache")
    _write_metadata(generator, sdf_path, "BBBBB", {"A": "SBMA", "B": "PEGMA"})
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.0, 1.0],
        length=5,
        type_prefix="PEGMA",
        sdf_directory=tmp_path,
        cache_directory=tmp_path / "cache",
        generation_mode="dynamic",
        monomer_smiles={"SBMA": "C", "PEGMA": "CO"},
        monomer_names={"A": "SBMA", "B": "PEGMA"},
        reactions=object(),
    )
    builder._fragment_generator = object()
    builder._polymer_generator = generator

    monkeypatch.setattr(builder, "_load_from_sdf", lambda path: sentinel)
    monkeypatch.setattr(
        builder,
        "_generate_polymer",
        lambda sequence: pytest.fail("Validated generated SDF was not loaded"),
    )

    assert builder._get_or_create_molecule("BBBBB") is sentinel


def test_dynamic_builder_regenerates_generated_sdf_directory_hit_with_stale_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should regenerate generated-style SDF hits with stale metadata."""
    sdf_path = tmp_path / "PEGMA_seq=BBBBB_5-mer_charged.sdf"
    sdf_path.write_text("stale generated artifact", encoding="utf-8")
    current_group = FakeMonomerGroup(
        {
            "SBMA_1-site": ["sbma-terminal"],
            "SBMA_2-site": ["sbma-middle"],
            "PEGMA_1-site": ["pegma-terminal"],
            "PEGMA_2-site": ["pegma-middle"],
        }
    )
    stale_group = FakeMonomerGroup(
        {
            "SBMA_1-site": ["sbma-terminal"],
            "SBMA_2-site": ["sbma-middle"],
            "PEGMA_1-site": ["pegma-terminal"],
            "PEGMA_2-site": ["stale-pegma-middle"],
        }
    )
    _write_metadata(
        PolymerGenerator(monomer_group=stale_group, cache_directory=tmp_path / "stale"),
        sdf_path,
        "BBBBB",
        {"A": "SBMA", "B": "PEGMA"},
    )
    sentinel = object()
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.0, 1.0],
        length=5,
        type_prefix="PEGMA",
        sdf_directory=tmp_path,
        cache_directory=tmp_path / "cache",
        generation_mode="dynamic",
        monomer_smiles={"SBMA": "C", "PEGMA": "CO"},
        monomer_names={"A": "SBMA", "B": "PEGMA"},
        reactions=object(),
    )
    builder._fragment_generator = object()
    builder._polymer_generator = PolymerGenerator(
        monomer_group=current_group,
        cache_directory=tmp_path / "cache",
    )

    monkeypatch.setattr(
        builder,
        "_load_from_sdf",
        lambda path: pytest.fail("Stale generated SDF was loaded directly"),
    )
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("BBBBB") is sentinel


def test_dynamic_filename_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Filename construction should reject unknown labels with ValueError."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator._make_polymer_filename("AXA", {"A": "EGMA"})


def test_dynamic_structure_build_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Structure building should reject unknown labels before direct mapping lookup."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator._build_polymer_structure("AXA", {"A": "EGMA"})


def test_dynamic_generation_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Generation should reject unknown labels before importing OpenFF topology."""
    generator = _make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator.generate_polymer("AXA", {"A": "EGMA"})


def test_cached_builder_loads_sdf_directory_hit_without_dynamic_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Cached mode should continue to load user-provided SDF directory hits."""
    sdf_path = tmp_path / "EGMA_seq=A_1-mer_charged.sdf"
    sdf_path.write_text("curated user structure", encoding="utf-8")
    sentinel = object()
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=1,
        type_prefix="EGMA",
        sdf_directory=tmp_path,
        cache_directory=tmp_path / "cache",
        generation_mode="cached",
    )

    monkeypatch.setattr(builder, "_load_from_sdf", lambda path: sentinel)

    assert builder._get_or_create_molecule("A") is sentinel
