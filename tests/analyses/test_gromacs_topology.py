"""Reading GROMACS runs whose prod.tpr MDAnalysis cannot parse.

The fixture in ``tests/data/gromacs/tpr_fallback`` is two methanols and five
SETTLE waters split into two blocks by a sodium ion, compiled by grompp from
GROMACS 2025.4 (``system_gmx2025.tpr``, tpx 137) and 2026.0
(``system_gmx2026.tpr``, tpx 138). Its atom types give atomic numbers in
both ``[ atomtypes ]`` layouts, with and without a bonded type.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")

from polyzymd.analyses.exceptions import ProtocolError  # noqa: E402
from polyzymd.analyses.shared import gromacs, loader  # noqa: E402
from polyzymd.analyses.shared.topology import topology_bond_source  # noqa: E402

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "gromacs" / "tpr_fallback"
ATOM_ATTRIBUTES = (
    "ids",
    "names",
    "types",
    "resids",
    "resnums",
    "resnames",
    "segids",
    "chainIDs",
    "moltypes",
    "molnums",
    "masses",
    "charges",
    "elements",
)


@pytest.fixture
def run_dir(tmp_path: Path) -> Path:
    """A copy of the fixture laid out as a PolyzyMD GROMACS run, with the 2026 TPR as prod.tpr."""
    folder = tmp_path / "run_1" / "gromacs"
    folder.mkdir(parents=True)
    for name in ("system.top", "system.gro") + tuple(p.name for p in FIXTURE.glob("*.itp")):
        shutil.copy(FIXTURE / name, folder / name)
    shutil.copy(FIXTURE / "system_gmx2026.tpr", folder / "prod.tpr")
    return folder


def _reads(path: Path) -> bool:
    try:
        mda.Universe(str(path))
    except ValueError:
        return False
    return True


def test_atomic_numbers_read_both_atomtypes_layouts() -> None:
    assert gromacs.atomic_numbers(FIXTURE / "system.top") == {
        "AT_0": 6,
        "AT_1": 1,
        "AT_2": 8,
        "AT_3": 1,
        "AT_4": 8,
        "AT_5": 1,
        "AT_6": 11,
    }


def test_atomic_numbers_skip_types_without_one(tmp_path: Path) -> None:
    top = tmp_path / "x.top"
    top.write_text(
        '#include "missing_forcefield.itp"\n[ atomtypes ]\nCT  12.011  0.0  A  0.3  0.4\n'
        "HC  HX  1.008  0.0  A  0.2  0.1\nOW  8  15.999  0.0  A  0.3  0.6\n"
    )
    assert gromacs.atomic_numbers(top) == {"OW": 8}


def test_top_universe_is_laid_out_as_mdanalysis_reads_the_tpr() -> None:
    tpr = FIXTURE / "system_gmx2025.tpr"
    if not _reads(tpr):
        pytest.skip("this MDAnalysis cannot read GROMACS 2025 TPR files")
    expected = mda.Universe(str(tpr))
    built = gromacs.universe_from_gromacs_top(FIXTURE / "system.top", FIXTURE / "system.gro")
    for name in ATOM_ATTRIBUTES:
        want, got = (
            np.asarray(getattr(expected.atoms, name)),
            np.asarray(getattr(built.atoms, name)),
        )
        assert got.dtype == want.dtype, name
        np.testing.assert_array_equal(got, want, err_msg=name)
    np.testing.assert_array_equal(built.segments.segids, expected.segments.segids)
    for group in ("bonds", "angles", "dihedrals", "impropers"):
        want = {tuple(i) for i in getattr(expected, group).indices}
        assert {tuple(i) for i in getattr(built, group).indices} == want, group


def test_top_universe_layout() -> None:
    u = gromacs.universe_from_gromacs_top(FIXTURE / "system.top", FIXTURE / "system.gro")
    assert list(u.segments.segids) == ["seg_0_MOL0", "seg_1_SOL", "seg_2_MOL1", "seg_3_SOL"]
    assert list(u.residues.resids) == list(range(1, 9))
    assert set(u.select_atoms("resname HOH").elements) == {"O", "H"}
    # Each SETTLE water has its two O-H bonds; methanol keeps its five bonds.
    water = u.select_atoms("resname HOH")
    assert len(water.bonds) == 2 * len(water.residues)
    assert len(u.select_atoms("resname MOH").bonds) == 10
    assert u.atoms.charges.dtype == np.float64
    assert len(u.trajectory) == 1


def test_tpr_unsupported_recognises_only_the_version_refusal() -> None:
    try:
        raise ValueError("Failed") from None
    except ValueError as plain:
        assert not gromacs.tpr_unsupported(plain)
    try:
        try:
            raise NotImplementedError("Your tpx version is 138")
        except NotImplementedError:
            raise ValueError("Failed to construct topology")
    except ValueError as wrapped:
        assert gromacs.tpr_unsupported(wrapped)


def test_open_universe_reads_the_top_when_the_tpr_is_too_new(run_dir: Path, caplog) -> None:
    tpr = run_dir / "prod.tpr"
    if _reads(tpr):
        pytest.skip("this MDAnalysis reads GROMACS 2026 TPR files")
    loader._WARNED_TPR_FALLBACK_PATHS.discard(tpr)
    with caplog.at_level("WARNING"):
        u = loader.open_universe(tpr, [run_dir / "system.gro"])
    assert "cannot read" in caplog.text and "system.top" in caplog.text
    assert u._polyzymd_bond_source == "top"
    assert topology_bond_source(u) == (True, "top")
    assert len(u.atoms) == 28
    assert list(u.segments.segids)[0] == "seg_0_MOL0"


def test_open_universe_needs_one_top(run_dir: Path) -> None:
    tpr = run_dir / "prod.tpr"
    if _reads(tpr):
        pytest.skip("this MDAnalysis reads GROMACS 2026 TPR files")
    (run_dir / "system.top").unlink()
    with pytest.raises(ProtocolError, match="exactly one"):
        loader.open_universe(tpr, [run_dir / "system.gro"])


def test_open_universe_reports_tpr_bonds(tmp_path: Path) -> None:
    tpr = FIXTURE / "system_gmx2025.tpr"
    if not _reads(tpr):
        pytest.skip("this MDAnalysis cannot read GROMACS 2025 TPR files")
    u = loader.open_universe(tpr, [FIXTURE / "system.gro"])
    assert topology_bond_source(u) == (True, "tpr")


def test_open_universe_propagates_other_errors(tmp_path: Path) -> None:
    broken = tmp_path / "prod.tpr"
    broken.write_bytes(b"not a tpr")
    with pytest.raises(Exception) as caught:
        loader.open_universe(broken, [])
    assert not isinstance(caught.value, ProtocolError)


def _write_build_pdb(path: Path, chains: list[str]) -> None:
    u = gromacs.universe_from_gromacs_top(FIXTURE / "system.top", FIXTURE / "system.gro")
    u.atoms.chainIDs = np.asarray(chains, dtype=object)
    u.atoms.write(str(path))


def test_build_chain_ids_replace_moltype_chains(run_dir: Path) -> None:
    chains = ["C"] * 12 + ["D"] * 16
    _write_build_pdb(run_dir.parent / "solvated_system.pdb", chains)
    assert gromacs.build_pdb_file(run_dir / "prod.tpr") == run_dir.parent / "solvated_system.pdb"
    u = loader.open_universe(run_dir / "prod.tpr", [run_dir / "system.gro"])
    assert u._polyzymd_chain_ids["applied"]
    assert list(u.atoms.chainIDs) == chains
    assert len(u.select_atoms("chainID C")) == 12


def test_build_chain_ids_need_matching_residues(tmp_path: Path) -> None:
    u = gromacs.universe_from_gromacs_top(FIXTURE / "system.top", FIXTURE / "system.gro")
    pdb = tmp_path / "solvated_system.pdb"
    u.select_atoms("not resname Na+").write(str(pdb))
    metadata = gromacs.apply_build_chain_ids(u, pdb)
    assert not metadata["applied"] and "atoms" in metadata["reason"]
    assert u.atoms.chainIDs[0] == "MOL0"
    assert gromacs.apply_build_chain_ids(u, None)["reason"] == "no solvated_system.pdb"
