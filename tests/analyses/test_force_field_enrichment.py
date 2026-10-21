"""Tests for openmm_system_file, read_openmm_system and enrich_universe_force_field.

The systems are serialized by OpenMM itself (tests/_support/openmm_system.py),
so the reader sees the layout OpenMM writes. The study tests write an OpenMM
run directory (tests/_support/analysis_testkit.py) of a serine side chain,
an SBM atom and a rigid water, whose PDB records one wrong bond, and put the
system XML beside its trajectory; the loaded universe must carry the
system's charges and bonds instead.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import polyzymd as pz
from polyzymd.analyses.shared.loader import (
    enrich_universe_force_field,
    openmm_system_file,
    read_openmm_system,
)
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config
from tests._support.openmm_system import openmm_system_xml, write_openmm_system

mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("openmm")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

#: SER 12 OG, HG and CB; SBM 1 O1; HOH 100 OW, HW1 and HW2.
NAMES = ["OG", "HG", "CB", "O1", "OW", "HW1", "HW2"]
ELEMENTS = ["O", "H", "C", "O", "O", "H", "H"]
RESINDEX = [0, 0, 0, 1, 2, 2, 2]
RESNAMES = ["SER", "SBM", "HOH"]
RESIDS = [12, 1, 100]
CHAINS = ["A", "A", "A", "C", "W", "W", "W"]
CHARGES = [-0.65, 0.42, 0.23, -0.5, -0.834, 0.417, 0.417]
HARMONIC = [(0, 2)]
#: OpenMM holds bonds to hydrogen and rigid water, H-H included, as constraints.
CONSTRAINTS = [(0, 1), (4, 5), (4, 6), (5, 6)]
FORCE_FIELD_BONDS = {(0, 1), (0, 2), (4, 5), (4, 6)}
#: A bond the PDB records and the system does not have.
PDB_BONDS = [(2, 3)]


def _bond_set(universe) -> set[tuple[int, int]]:
    return {tuple(sorted(int(i) for i in bond.indices)) for bond in universe.bonds}


def _universe(bonds=None, elements=ELEMENTS) -> "mda.Universe":
    universe = mda.Universe.empty(
        len(NAMES), n_residues=len(RESNAMES), atom_resindex=RESINDEX, trajectory=True
    )
    universe.add_TopologyAttr("names", NAMES)
    universe.add_TopologyAttr("elements", list(elements))
    universe.add_TopologyAttr("resnames", RESNAMES)
    universe.add_TopologyAttr("resids", RESIDS)
    if bonds is not None:
        universe.add_TopologyAttr("bonds", bonds)
    return universe


# ---------------------------------------------------------------------------
# openmm_system_file
# ---------------------------------------------------------------------------


def test_openmm_system_file_finds_the_system_beside_a_segment_trajectory(tmp_path) -> None:
    segment = tmp_path / "production_3"
    segment.mkdir()
    trajectory = segment / "production_3_trajectory.dcd"
    trajectory.touch()
    system = segment / "production_3_system.xml"

    assert openmm_system_file(trajectory) is None
    system.write_text("<System/>")
    assert openmm_system_file(trajectory) == system
    assert openmm_system_file(str(trajectory)) == system


def test_openmm_system_file_uses_the_whole_stem_without_the_trajectory_suffix(tmp_path) -> None:
    (tmp_path / "run_system.xml").write_text("<System/>")
    (tmp_path / "production_0_system.xml").write_text("<System/>")

    assert openmm_system_file(tmp_path / "run.dcd") == tmp_path / "run_system.xml"
    assert openmm_system_file(tmp_path / "other_trajectory.dcd") is None


def test_openmm_system_file_ignores_a_directory_of_that_name(tmp_path) -> None:
    (tmp_path / "production_0_system.xml").mkdir()

    assert openmm_system_file(tmp_path / "production_0_trajectory.dcd") is None


# ---------------------------------------------------------------------------
# read_openmm_system
# ---------------------------------------------------------------------------


def test_read_openmm_system_returns_nonbonded_charges_and_bonds_with_constraints(
    tmp_path,
) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES, HARMONIC, CONSTRAINTS))

    charges, bonds = read_openmm_system(path)

    assert charges.dtype == np.float64
    assert charges.tolist() == pytest.approx(CHARGES)
    assert sorted(bonds) == sorted(HARMONIC + CONSTRAINTS)


def test_read_openmm_system_ignores_other_forces_with_particles_or_bonds(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES, HARMONIC, CONSTRAINTS, decoys=True))

    charges, bonds = read_openmm_system(path)

    # The GBSAOBCForce charges of 9 and the CustomBondForce bond (0, 6) are left out.
    assert charges.tolist() == pytest.approx(CHARGES)
    assert sorted(bonds) == sorted(HARMONIC + CONSTRAINTS)


def test_read_openmm_system_of_a_system_without_bonds(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml([0.5, -0.5]))

    charges, bonds = read_openmm_system(path)

    assert charges.tolist() == [0.5, -0.5]
    assert bonds == []


# ---------------------------------------------------------------------------
# enrich_universe_force_field
# ---------------------------------------------------------------------------


def test_enrichment_adds_charges_and_replaces_bonds_without_h_h_constraints(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES, HARMONIC, CONSTRAINTS))
    universe = _universe(bonds=PDB_BONDS)

    metadata = enrich_universe_force_field(universe, path)

    assert universe.atoms.charges.tolist() == pytest.approx(CHARGES)
    assert _bond_set(universe) == FORCE_FIELD_BONDS
    assert [a.name for a in universe.atoms[5].bonded_atoms] == ["OW"]
    assert metadata == {"applied": True, "source": str(path), "bonds": len(FORCE_FIELD_BONDS)}
    assert universe._polyzymd_force_field is metadata


def test_enrichment_adds_bonds_to_a_universe_that_had_none(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES, HARMONIC, CONSTRAINTS))
    universe = _universe()

    enrich_universe_force_field(universe, path)

    assert _bond_set(universe) == FORCE_FIELD_BONDS


def test_enrichment_accepts_lowercase_or_padded_elements(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES, HARMONIC, CONSTRAINTS))
    universe = _universe(elements=[" o", "h ", "c", "O", "o", " H", "h"])

    enrich_universe_force_field(universe, path)

    assert _bond_set(universe) == FORCE_FIELD_BONDS


def test_enrichment_without_a_system_file_leaves_the_universe_unchanged() -> None:
    universe = _universe(bonds=PDB_BONDS)

    metadata = enrich_universe_force_field(universe, None)

    assert metadata == {"applied": False, "source": None, "reason": "no system XML"}
    assert universe._polyzymd_force_field is metadata
    assert _bond_set(universe) == set(PDB_BONDS)
    assert not hasattr(universe.atoms, "charges")


def test_enrichment_with_another_particle_count_leaves_the_universe_unchanged(tmp_path) -> None:
    path = tmp_path / "system.xml"
    path.write_text(openmm_system_xml(CHARGES[:4], HARMONIC, CONSTRAINTS[:1]))
    universe = _universe(bonds=PDB_BONDS)

    metadata = enrich_universe_force_field(universe, path)

    assert metadata == {"applied": False, "source": str(path), "reason": "4 particles for 7 atoms"}
    assert universe._polyzymd_force_field is metadata
    assert _bond_set(universe) == set(PDB_BONDS)
    assert not hasattr(universe.atoms, "charges")


# ---------------------------------------------------------------------------
# Through the study loader
# ---------------------------------------------------------------------------

COORDINATES = np.array(
    [
        [
            [10.0, 10.0, 10.0],
            [10.97, 10.0, 10.0],
            [9.5, 11.4, 10.0],
            [13.9, 10.0, 10.0],
            [30.0, 30.0, 30.0],
            [30.96, 30.0, 30.0],
            [29.76, 30.93, 30.0],
        ]
    ]
    * 3,
    dtype=np.float32,
)


def _config(tmp_path: Path, *, system: bool = True, charges=CHARGES) -> Path:
    config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "A" / "scratch")
    for replicate in (1, 2):
        run_dir = write_openmm_frames(
            config,
            replicate,
            COORDINATES,
            RESINDEX,
            resids=RESIDS,
            names=NAMES,
            resnames=RESNAMES,
            elements=ELEMENTS,
            chain_ids=CHAINS,
            dimensions=[40.0, 40.0, 40.0, 90.0, 90.0, 90.0],
            bonds=PDB_BONDS,
        )
        if system:
            write_openmm_system(run_dir, charges, HARMONIC, CONSTRAINTS)
    return config


def test_the_study_universe_carries_the_charges_and_bonds_of_the_run_system(tmp_path) -> None:
    study = pz.Study.from_configs({"A": _config(tmp_path)}, equilibration="0ns")

    for replicate in study["A"].replicates:
        universe = replicate.universe()
        metadata = universe._polyzymd_force_field
        assert metadata["applied"] is True
        assert Path(metadata["source"]).name == "production_0_system.xml"
        assert Path(metadata["source"]).parent.parent.name.endswith(f"_run{replicate.index}")
        assert metadata["bonds"] == len(FORCE_FIELD_BONDS)
        assert universe.atoms.charges.tolist() == pytest.approx(CHARGES, abs=1e-6)
        assert _bond_set(universe) == FORCE_FIELD_BONDS
        assert universe.select_atoms("resname HOH and element H").n_atoms == 2


def test_the_study_universe_keeps_the_pdb_bonds_without_a_system_file(tmp_path) -> None:
    study = pz.Study.from_configs({"A": _config(tmp_path, system=False)}, equilibration="0ns")

    universe = study["A"].replicates[0].universe()

    assert universe._polyzymd_force_field == {
        "applied": False,
        "source": None,
        "reason": "no system XML",
    }
    assert _bond_set(universe) == set(PDB_BONDS)


def test_the_study_universe_ignores_a_system_of_another_particle_count(tmp_path) -> None:
    config = _config(tmp_path, charges=CHARGES + [0.0])
    study = pz.Study.from_configs({"A": config}, equilibration="0ns")

    universe = study["A"].replicates[0].universe()

    assert universe._polyzymd_force_field["applied"] is False
    assert universe._polyzymd_force_field["reason"] == "8 particles for 7 atoms"
    assert _bond_set(universe) == set(PDB_BONDS)
