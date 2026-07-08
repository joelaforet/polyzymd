"""Pure PolyzyMD unit tests for dynamic polymer generation mapping."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from _polymer_generator_helpers import (
    FakeMonomerGroup,
    make_generator,
    make_list_fragment_group,
    patch_structure_build,
    write_metadata,
)

import polyzymd.builders.polymer_generator as polymer_generator_module
from polyzymd.builders.polymer_generator import PolymerGenerator


def test_two_monomer_b_sequence_uses_explicit_pegma_middle_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """BBBBB should pass BBB and map B to the PEGMA 2-site fragment."""
    calls: list[dict[str, object]] = []
    patch_structure_build(monkeypatch, calls)
    generator = make_generator(
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
    patch_structure_build(monkeypatch, calls)
    generator = make_generator(
        tmp_path,
        ["SBMA_1-site", "SBMA_2-site", "PEGMA_1-site", "PEGMA_2-site"],
    )

    generator._build_polymer_structure("BBBBB", {"A": "SBMA", "B": "PEGMA"})

    monomer_group = calls[0]["monomers"]
    assert isinstance(monomer_group, FakeMonomerGroup)
    assert monomer_group.term_orient == {"head": "PEGMA_1-site", "tail": "PEGMA_1-site"}


def test_single_monomer_sequence_uses_explicit_egma_middle_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """AAAAA should pass an explicit EGMA mapping instead of Polymerist defaults."""
    calls: list[dict[str, object]] = []
    patch_structure_build(monkeypatch, calls)
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    generator._build_polymer_structure("AAAAA", {"A": "EGMA"})

    assert calls[0]["sequence"] == "AAA"
    assert calls[0]["sequence_map"] == {"A": "EGMA_2-site"}


def test_stale_dynamic_cache_without_metadata_is_not_loaded(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A dynamic charged SDF without metadata should not be directly reused."""
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")

    def fail_load_cache(path: Path) -> object:
        """Fail if stale SDF loading is attempted.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
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
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    generator._get_dynamic_cache_metadata_path(cache_path).write_text(
        json.dumps(metadata),
        encoding="utf-8",
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if non-object metadata cache loading is attempted.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
        pytest.fail(f"Non-object metadata dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_dynamic_cache_rejects_charger_mismatch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache metadata must include and validate the charger type."""
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    write_metadata(
        generator,
        cache_path,
        "AAAAA",
        {"A": "EGMA"},
        overrides={"charger_type": "espaloma"},
    )

    def fail_load_cache(path: Path) -> object:
        """Fail if mismatched dynamic cache loading is attempted.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
        pytest.fail(f"Mismatched dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_dynamic_cache_rejects_fragment_fingerprint_mismatch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic cache metadata must validate MonomerGroup fragment chemistry."""
    current_generator = PolymerGenerator(
        monomer_group=make_list_fragment_group(),
        cache_directory=tmp_path,
    )
    stale_generator = PolymerGenerator(
        monomer_group=make_list_fragment_group({"EGMA_2-site": ["[C:1]-[C:2]", "[N:3]"]}),
        cache_directory=tmp_path,
    )
    cache_path = tmp_path / "EGMA_seq=AAAAA_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    write_metadata(stale_generator, cache_path, "AAAAA", {"A": "EGMA"})

    def fail_load_cache(path: Path) -> object:
        """Fail if fingerprint-mismatched cache loading is attempted.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
        pytest.fail(f"Fingerprint-mismatched dynamic cache was loaded directly: {path}")

    monkeypatch.setattr(polymer_generator_module, "_topology_from_sdf", fail_load_cache)

    assert current_generator.get_cached_polymer("AAAAA", {"A": "EGMA"}) is None


def test_list_fragment_contents_affect_monomer_group_fingerprint(tmp_path: Path) -> None:
    """Same fragment keys with different list contents should not share a digest."""
    first_generator = PolymerGenerator(
        monomer_group=make_list_fragment_group(),
        cache_directory=tmp_path / "first",
    )
    second_generator = PolymerGenerator(
        monomer_group=make_list_fragment_group({"EGMA_2-site": ["[C:1]-[C:2]", "[N:3]"]}),
        cache_directory=tmp_path / "second",
    )

    first_fingerprint = first_generator._build_monomer_group_fingerprint()
    second_fingerprint = second_generator._build_monomer_group_fingerprint()

    assert first_fingerprint["algorithm"] == "polyzymd-monomer-group-sha256-v2"
    assert first_fingerprint["digest"] != second_fingerprint["digest"]


def test_dynamic_generator_rejects_short_sequences(tmp_path: Path) -> None:
    """Lower-level dynamic generation should reject sequences shorter than three."""
    generator = PolymerGenerator(monomer_group=make_list_fragment_group(), cache_directory=tmp_path)

    with pytest.raises(ValueError, match="requires sequence length >= 3"):
        generator.generate_polymer("AA", {"A": "EGMA"})


def test_dynamic_filename_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Filename construction should reject unknown labels with ValueError."""
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator._make_polymer_filename("AXA", {"A": "EGMA"})


def test_dynamic_structure_build_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Structure building should reject unknown labels before direct mapping lookup."""
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator._build_polymer_structure("AXA", {"A": "EGMA"})


def test_dynamic_generation_rejects_unknown_sequence_label(tmp_path: Path) -> None:
    """Generation should reject unknown labels before importing OpenFF topology."""
    generator = make_generator(tmp_path, ["EGMA_1-site", "EGMA_2-site"])

    with pytest.raises(ValueError, match="No monomer name configured for sequence label: X"):
        generator.generate_polymer("AXA", {"A": "EGMA"})
