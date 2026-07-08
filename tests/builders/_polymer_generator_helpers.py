"""Shared helpers for polymer generator tests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import polyzymd.builders.polymer_generator as polymer_generator_module
from polyzymd.builders.polymer_generator import PolymerGenerator


class FakeMonomerGroup:
    """Minimal MonomerGroup stand-in for local unit tests."""

    def __init__(self, monomers: dict[str, object]) -> None:
        """Initialize the fake group with named fragments.

        Parameters
        ----------
        monomers : dict[str, object]
            Fragment payloads keyed by Polymerist-style fragment names.
        """
        self.monomers = monomers
        self.term_orient: dict[str, str] = {}


def make_generator(tmp_path: Path, fragment_names: list[str]) -> PolymerGenerator:
    """Create a PolymerGenerator with fake fragment objects.

    Parameters
    ----------
    tmp_path : Path
        Temporary cache directory.
    fragment_names : list[str]
        Fragment names to expose through the fake MonomerGroup.

    Returns
    -------
    PolymerGenerator
        Generator configured with a fake MonomerGroup.
    """
    monomer_group = FakeMonomerGroup({name: [f"{name}:smarts"] for name in fragment_names})
    return PolymerGenerator(monomer_group=monomer_group, cache_directory=tmp_path)


def make_list_fragment_group(overrides: dict[str, list[str]] | None = None) -> FakeMonomerGroup:
    """Create a fake MonomerGroup with Polymerist-like list fragment entries.

    Parameters
    ----------
    overrides : dict[str, list[str]] | None, optional
        Fragment payloads to merge into the default fake fragments.

    Returns
    -------
    FakeMonomerGroup
        Fake group with deterministic list-backed fragment values.
    """
    fragments = {
        "EGMA_1-site": ["[C:1]=[C:2]", "[O:3]"],
        "EGMA_2-site": ["[C:1]-[C:2]", "[O:3]"],
    }
    fragments.update(overrides or {})
    return FakeMonomerGroup(fragments)


def patch_structure_build(
    monkeypatch: pytest.MonkeyPatch,
    calls: list[dict[str, object]],
) -> None:
    """Patch Polymerist structure helpers to capture PolyzyMD boundary calls.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest monkeypatch fixture.
    calls : list[dict[str, object]]
        Mutable list that receives keyword arguments passed to the build wrapper.
    """

    def fake_build_linear_polymer(**kwargs: Any) -> object:
        """Capture the Polymerist build arguments.

        Parameters
        ----------
        **kwargs : Any
            Keyword arguments from PolyzyMD's Polymerist wrapper.

        Returns
        -------
        object
            Sentinel chain object.
        """
        calls.append(kwargs)
        return object()

    def fake_write_pdb(pdb_path: Path, chain: object, resname_map: dict[str, str]) -> None:
        """Create the expected PDB output without invoking OpenMM writers.

        Parameters
        ----------
        pdb_path : Path
            Destination PDB path.
        chain : object
            Ignored chain sentinel.
        resname_map : dict[str, str]
            Ignored residue-name mapping.
        """
        pdb_path.write_text("MODEL\nEND\n", encoding="utf-8")

    monkeypatch.setattr(polymer_generator_module, "_make_monomer_group", FakeMonomerGroup)
    monkeypatch.setattr(
        polymer_generator_module, "_build_linear_polymer", fake_build_linear_polymer
    )
    monkeypatch.setattr(polymer_generator_module, "_mbmol_to_rdmol", lambda chain: object())
    monkeypatch.setattr(polymer_generator_module, "_summarize_ring_piercing", lambda mol: {})
    monkeypatch.setattr(polymer_generator_module, "_mbmol_to_openmm_pdb", fake_write_pdb)


def write_metadata(
    generator: PolymerGenerator,
    sdf_path: Path,
    sequence: str,
    monomer_names: dict[str, str],
    overrides: dict[str, object] | None = None,
) -> None:
    """Write dynamic cache metadata with optional top-level overrides.

    Parameters
    ----------
    generator : PolymerGenerator
        Generator used to build the expected metadata payload.
    sdf_path : Path
        Charged SDF path whose sidecar should be written.
    sequence : str
        Polymer sequence for metadata validation.
    monomer_names : dict[str, str]
        Mapping from sequence labels to monomer names.
    overrides : dict[str, object] | None, optional
        Values to merge into the generated metadata payload.
    """
    metadata = generator._build_dynamic_cache_metadata(sequence, monomer_names)
    metadata.update(overrides or {})
    generator._get_dynamic_cache_metadata_path(sdf_path).write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
