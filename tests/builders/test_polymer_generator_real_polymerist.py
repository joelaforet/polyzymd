"""Regression tests that exercise real Polymerist/RDKit polymer data paths."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import polyzymd.builders.polymer_generator as polymer_generator_module
from polyzymd.builders.polymer_generator import PolymerGenerator

pytest.importorskip("polymerist")
pytest.importorskip("rdkit")


def _real_halogen_fragments() -> dict[str, list[Any]]:
    """Create valid Polymerist SMARTS fragments with PolyzyMD monomer names.

    Returns
    -------
    dict[str, list[Any]]
        Polymerist fragment payloads keyed by PolyzyMD fragment names.
    """
    from polymerist.polymers.monomers.fragments import HALOGENATED_HYDROCARBON_FRAGMENTS

    return {
        "SBMA_1-site": HALOGENATED_HYDROCARBON_FRAGMENTS["chlor_term_1"],
        "SBMA_2-site": HALOGENATED_HYDROCARBON_FRAGMENTS["chlor_mid_1"],
        "PEGMA_1-site": HALOGENATED_HYDROCARBON_FRAGMENTS["brom_term_1"],
        "PEGMA_2-site": HALOGENATED_HYDROCARBON_FRAGMENTS["brom_mid_1"],
    }


def _make_real_halogen_monomer_group(overrides: dict[str, list[Any]] | None = None) -> Any:
    """Create a real Polymerist MonomerGroup backed by list SMARTS data.

    Parameters
    ----------
    overrides : dict[str, list[Any]] | None, optional
        Fragment entries to merge into the default real fragments.

    Returns
    -------
    Any
        Polymerist MonomerGroup instance.
    """
    from polymerist.polymers.monomers import MonomerGroup

    fragments = _real_halogen_fragments()
    fragments.update(overrides or {})
    return MonomerGroup(monomers=fragments)


def _count_chain_halogens(chain: object) -> dict[str, int]:
    """Count diagnostic halogens in a real Polymerist-built chain.

    Parameters
    ----------
    chain : object
        Polymerist chain object to inspect.

    Returns
    -------
    dict[str, int]
        Counts of bromine and chlorine atoms.
    """
    from polymerist.polymers.building import mbmol_to_rdmol

    mol = mbmol_to_rdmol(chain)
    return {
        "Br": sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Br"),
        "Cl": sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl"),
    }


def _write_metadata(
    generator: PolymerGenerator,
    sdf_path: Path,
    sequence: str,
    monomer_names: dict[str, str],
) -> None:
    """Write dynamic cache metadata for real Polymerist tests.

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
    """
    metadata = generator._build_dynamic_cache_metadata(sequence, monomer_names)
    generator._get_dynamic_cache_metadata_path(sdf_path).write_text(
        polymer_generator_module.json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_real_polymerist_sequence_map_keeps_b_middle_on_pegma_fragment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """BBBBB should build with PEGMA middle fragments through real Polymerist."""
    from polymerist.polymers.building import build_linear_polymer

    def fake_write_pdb(pdb_path: Path, chain: object, resname_map: dict[str, str]) -> None:
        """Create the expected PDB output after the real chain build.

        Parameters
        ----------
        pdb_path : Path
            Destination PDB path.
        chain : object
            Polymerist chain object.
        resname_map : dict[str, str]
            Residue-name mapping passed by PolyzyMD.
        """
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

    chain, pdb_path = generator._build_polymer_structure("BBBBB", {"A": "SBMA", "B": "PEGMA"})

    assert pdb_path.exists()
    assert _count_chain_halogens(chain) == {"Br": 5, "Cl": 0}


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

    assert current_generator.get_cached_polymer("BBBBB", {"A": "SBMA", "B": "PEGMA"}) is None
