"""Tests for PolymerBuilder dynamic cache routing."""

from __future__ import annotations

from pathlib import Path

import pytest
from _polymer_generator_helpers import FakeMonomerGroup, write_metadata

from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.polymer_generator import PolymerGenerator


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
        """Fail if PolymerBuilder bypasses dynamic cache validation.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
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
        """Fail if generated-style SDF loading bypasses validation.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
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
    write_metadata(generator, sdf_path, "BBBBB", {"A": "SBMA", "B": "PEGMA"})
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
    write_metadata(
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
