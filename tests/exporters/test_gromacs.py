"""Tests for GROMACS exporter bug fixes.

Covers:
- Bug 1: Position restraints applied to ALL polymer ITP files
  (not just the first one)
- Bug 2: GRO residue numbering for multi-residue molecule copies
  (globally sequential instead of sliding offset)
- ITP parsing helpers (_get_heavy_atom_indices_from_itp,
  _get_atom_count_from_itp, _get_molecule_name_from_itp,
  _extract_mol_index, _find_all_polymer_itps)
- GRO residue fix helpers (_parse_molecules_from_top,
  _parse_itp_atom_and_residue_counts, _fix_gro_residue_numbering)
"""

from importlib import resources
from pathlib import Path

import pytest

from polyzymd.exporters.gromacs import (
    GromacsExporter,
    GromacsRunner,
    MDPGenerator,
    PositionRestraintGenerator,
    RunScriptGenerator,
)


def _assert_no_unresolved_jinja(script: str) -> None:
    """Assert generated shell content has no unresolved Jinja syntax."""
    assert "{{" not in script
    assert "}}" not in script
    assert "{%" not in script
    assert "%}" not in script


# ---------------------------------------------------------------------------
# Fixtures — synthetic ITP content
# ---------------------------------------------------------------------------

# Minimal ITP for a polymer with 10 atoms (5 heavy + 5 hydrogen), 3 residues
POLYMER_ITP_TEMPLATE = """\
[ moleculetype ]
{mol_name}           3

[ atoms ]
;index, atom type, resnum, resname, name, cgnr, charge, mass
     1 {mol_name}_0       1 RES1     C1      1  -0.100000000000  12.010780000000
     2 {mol_name}_1       1 RES1     H1      2   0.050000000000   1.007940000000
     3 {mol_name}_2       1 RES1     O1      3  -0.200000000000  15.999430000000
     4 {mol_name}_3       2 RES2     C2      4  -0.100000000000  12.010780000000
     5 {mol_name}_4       2 RES2     H2      5   0.050000000000   1.007940000000
     6 {mol_name}_5       2 RES2     N1      6  -0.150000000000  14.006720000000
     7 {mol_name}_6       2 RES2     H3      7   0.050000000000   1.007940000000
     8 {mol_name}_7       3 RES3     C3      8  -0.100000000000  12.010780000000
     9 {mol_name}_8       3 RES3     H4      9   0.050000000000   1.007940000000
    10 {mol_name}_9       3 RES3     S1     10  -0.050000000000  32.065000000000

[ bonds ]
;ai    aj   funct r k
     1      2      1    0.10900  284512.000
"""

# Minimal ITP for protein with 6 atoms, 2 residues
PROTEIN_ITP_CONTENT = """\
[ moleculetype ]
MOL0           3

[ atoms ]
     1 MOL0_0       1 ALA      N       1  -0.300000000000  14.006720000000
     2 MOL0_1       1 ALA      CA      2   0.100000000000  12.010780000000
     3 MOL0_2       1 ALA      H       3   0.200000000000   1.007940000000
     4 MOL0_3       2 GLY      N       4  -0.300000000000  14.006720000000
     5 MOL0_4       2 GLY      CA      5   0.100000000000  12.010780000000
     6 MOL0_5       2 GLY      H       6   0.200000000000   1.007940000000

[ bonds ]
     1      2      1    0.14900  224262.400
"""

# Water ITP (3 atoms, 1 residue)
WATER_ITP_CONTENT = """\
[ moleculetype ]
water           2

[ atoms ]
     1 water_0       1 HOH      OW      1  -0.834000000000  15.999430000000
     2 water_1       1 HOH      HW1     2   0.417000000000   1.007940000000
     3 water_2       1 HOH      HW2     3   0.417000000000   1.007940000000

[ settles ]
     1      1    0.09572    0.15139
"""


def _write_itp(directory: Path, prefix: str, mol_name: str, content: str) -> Path:
    """Write an ITP file and return its path."""
    path = directory / f"{prefix}_{mol_name}.itp"
    path.write_text(content)
    return path


# ---------------------------------------------------------------------------
# Bug 1: ITP parsing helpers
# ---------------------------------------------------------------------------


class TestGetHeavyAtomIndicesFromItp:
    """_get_heavy_atom_indices_from_itp identifies non-H atoms by name."""

    def test_polymer_heavy_atoms(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL2", POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"))
        heavy = PositionRestraintGenerator._get_heavy_atom_indices_from_itp(itp)
        # Atoms: C1(1), H1(2), O1(3), C2(4), H2(5), N1(6), H3(7), C3(8), H4(9), S1(10)
        # Non-H: 1, 3, 4, 6, 8, 10
        assert heavy == [1, 3, 4, 6, 8, 10]

    def test_protein_heavy_atoms(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL0", PROTEIN_ITP_CONTENT)
        heavy = PositionRestraintGenerator._get_heavy_atom_indices_from_itp(itp)
        # N(1), CA(2), H(3), N(4), CA(5), H(6) -> non-H: 1, 2, 4, 5
        assert heavy == [1, 2, 4, 5]

    def test_water_heavy_atoms(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "water", WATER_ITP_CONTENT)
        heavy = PositionRestraintGenerator._get_heavy_atom_indices_from_itp(itp)
        # OW(1), HW1(2), HW2(3) -> non-H: 1 only
        # Wait — HW1 starts with 'H', so only OW is heavy
        assert heavy == [1]

    def test_nonexistent_file(self, tmp_path):
        path = tmp_path / "nope.itp"
        heavy = PositionRestraintGenerator._get_heavy_atom_indices_from_itp(path)
        assert heavy == []


class TestGetAtomCountFromItp:
    """_get_atom_count_from_itp counts all atoms in [ atoms ] section."""

    def test_polymer_atom_count(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL2", POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"))
        count = PositionRestraintGenerator._get_atom_count_from_itp(itp)
        assert count == 10

    def test_protein_atom_count(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL0", PROTEIN_ITP_CONTENT)
        count = PositionRestraintGenerator._get_atom_count_from_itp(itp)
        assert count == 6

    def test_water_atom_count(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "water", WATER_ITP_CONTENT)
        count = PositionRestraintGenerator._get_atom_count_from_itp(itp)
        assert count == 3


class TestGetMoleculeNameFromItp:
    """_get_molecule_name_from_itp extracts the moleculetype name."""

    def test_polymer_name(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL2", POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"))
        assert PositionRestraintGenerator._get_molecule_name_from_itp(itp) == "MOL2"

    def test_water_name(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "water", WATER_ITP_CONTENT)
        assert PositionRestraintGenerator._get_molecule_name_from_itp(itp) == "water"


class TestExtractMolIndex:
    """_extract_mol_index parses the MOL index from filenames."""

    def test_normal_cases(self):
        assert PositionRestraintGenerator._extract_mol_index(Path("sys_MOL0.itp")) == 0
        assert PositionRestraintGenerator._extract_mol_index(Path("sys_MOL12.itp")) == 12

    def test_no_match(self):
        assert PositionRestraintGenerator._extract_mol_index(Path("sys_water.itp")) is None


class TestFindAllPolymerItps:
    """_find_all_polymer_itps discovers polymer ITPs, skipping protein/ligand/water."""

    def _setup_system(self, tmp_path, has_protein=True, has_substrate=False):
        """Create a mock system with protein + 3 polymer types + water."""
        from unittest.mock import MagicMock

        prefix = "sys"

        # Create ITP files
        if has_protein:
            _write_itp(tmp_path, prefix, "MOL0", PROTEIN_ITP_CONTENT)

        mol_offset = (1 if has_protein else 0) + (1 if has_substrate else 0)
        for i in range(3):
            mol_name = f"MOL{mol_offset + i}"
            _write_itp(
                tmp_path,
                prefix,
                mol_name,
                POLYMER_ITP_TEMPLATE.format(mol_name=mol_name),
            )

        # Water ITP should be skipped (too few atoms)
        _write_itp(tmp_path, prefix, f"MOL{mol_offset + 3}", WATER_ITP_CONTENT)

        # Mock component_info
        component_info = MagicMock()
        component_info.has_protein = has_protein
        component_info.has_substrate = has_substrate

        gen = PositionRestraintGenerator.__new__(PositionRestraintGenerator)
        gen._component_info = component_info
        return gen, prefix

    def test_finds_all_polymer_itps_with_protein(self, tmp_path):
        gen, prefix = self._setup_system(tmp_path, has_protein=True)
        itps = gen._find_all_polymer_itps(tmp_path, prefix)
        names = [p.name for p in itps]
        assert names == ["sys_MOL1.itp", "sys_MOL2.itp", "sys_MOL3.itp"]

    def test_finds_all_polymer_itps_no_protein(self, tmp_path):
        gen, prefix = self._setup_system(tmp_path, has_protein=False)
        itps = gen._find_all_polymer_itps(tmp_path, prefix)
        names = [p.name for p in itps]
        assert names == ["sys_MOL0.itp", "sys_MOL1.itp", "sys_MOL2.itp"]

    def test_skips_water_itp(self, tmp_path):
        gen, prefix = self._setup_system(tmp_path, has_protein=True)
        itps = gen._find_all_polymer_itps(tmp_path, prefix)
        water_itps = [
            p
            for p in itps
            if "water" in (PositionRestraintGenerator._get_molecule_name_from_itp(p) or "").lower()
        ]
        assert len(water_itps) == 0

    def test_empty_directory(self, tmp_path):
        from unittest.mock import MagicMock

        component_info = MagicMock()
        component_info.has_protein = True
        component_info.has_substrate = False

        gen = PositionRestraintGenerator.__new__(PositionRestraintGenerator)
        gen._component_info = component_info
        itps = gen._find_all_polymer_itps(tmp_path, "sys")
        assert itps == []


class TestAppendPosresToItp:
    """_append_posres_to_itp correctly appends #ifdef POSRES block."""

    def test_appends_posres_block(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL2", POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"))
        PositionRestraintGenerator._append_posres_to_itp(
            itp, [1, 3, 4, 6, 8, 10], 1000.0, "POSRES_POLYMER", "polymer_heavy"
        )
        content = itp.read_text()

        assert "#ifdef POSRES_POLYMER" in content
        assert "[ position_restraints ]" in content
        assert "#endif" in content
        # Check all 6 heavy atom indices are present
        for idx in [1, 3, 4, 6, 8, 10]:
            assert f"{idx:6d}  1     1000.0    1000.0    1000.0" in content

    def test_does_not_duplicate_existing_content(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL2", POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"))
        original_content = itp.read_text()
        PositionRestraintGenerator._append_posres_to_itp(
            itp, [1], 500.0, "POSRES_POLYMER", "polymer_heavy"
        )
        content = itp.read_text()
        # Original content should still be there
        assert "[ moleculetype ]" in content
        assert "[ atoms ]" in content
        assert "[ bonds ]" in content


# ---------------------------------------------------------------------------
# Bug 2: GRO residue numbering helpers
# ---------------------------------------------------------------------------


class TestParseMoleculesFromTop:
    """_parse_molecules_from_top reads [ molecules ] section."""

    def test_normal_top_file(self, tmp_path):
        top = tmp_path / "sys.top"
        top.write_text("""\
; Topology file
#include "sys_MOL0.itp"
#include "sys_MOL1.itp"

[ system ]
test system

[ molecules ]
;name   number
MOL0    1
MOL1    4
MOL2    3
water   1000
NA      10
CL      10
""")
        result = GromacsExporter._parse_molecules_from_top(top)
        assert result == [
            ("MOL0", 1),
            ("MOL1", 4),
            ("MOL2", 3),
            ("water", 1000),
            ("NA", 10),
            ("CL", 10),
        ]

    def test_empty_molecules_section(self, tmp_path):
        top = tmp_path / "sys.top"
        top.write_text("[ system ]\ntest\n\n[ molecules ]\n")
        result = GromacsExporter._parse_molecules_from_top(top)
        assert result == []


class TestParseItpAtomAndResidueCounts:
    """_parse_itp_atom_and_residue_counts gets atom and residue counts."""

    def test_polymer_itp(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL1", POLYMER_ITP_TEMPLATE.format(mol_name="MOL1"))
        n_atoms, n_residues = GromacsExporter._parse_itp_atom_and_residue_counts(itp)
        assert n_atoms == 10
        assert n_residues == 3

    def test_protein_itp(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "MOL0", PROTEIN_ITP_CONTENT)
        n_atoms, n_residues = GromacsExporter._parse_itp_atom_and_residue_counts(itp)
        assert n_atoms == 6
        assert n_residues == 2

    def test_water_itp(self, tmp_path):
        itp = _write_itp(tmp_path, "sys", "water", WATER_ITP_CONTENT)
        n_atoms, n_residues = GromacsExporter._parse_itp_atom_and_residue_counts(itp)
        assert n_atoms == 3
        assert n_residues == 1


# ---------------------------------------------------------------------------
# Bug 2: Full GRO residue renumbering
# ---------------------------------------------------------------------------


def _make_gro_line(resid, resname, atomname, atomnum, x=0.0, y=0.0, z=0.0):
    """Build a single GRO atom line in the standard fixed-width format."""
    return f"{resid:5d}{resname:<5s}{atomname:>5s}{atomnum:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


class TestFixGroResidueNumbering:
    """_fix_gro_residue_numbering corrects multi-residue molecule copies."""

    def _setup_system(self, tmp_path):
        """Create a minimal system: 1 protein (2 res) + 2 polymer copies (3 res each) + 3 water.

        The broken GRO has the OpenFF copy_index offset bug:
        - Protein copy 0: resid 1, 1, 2, 2 (2 residues, correct)
        - Polymer copy 0: resid 1, 1, 2, 2, 3 (3 residues, correct)
        - Polymer copy 1: resid 2, 2, 3, 3, 4 (WRONG: offset by +1)
        - Water copy 0: resid 1 (correct, single residue)
        - Water copy 1: resid 2 (correct)
        - Water copy 2: resid 3 (correct)
        """
        prefix = "sys"

        # Write ITPs
        _write_itp(tmp_path, prefix, "MOL0", PROTEIN_ITP_CONTENT)
        _write_itp(
            tmp_path,
            prefix,
            "MOL1",
            POLYMER_ITP_TEMPLATE.format(mol_name="MOL1"),
        )
        # Water ITP with only 3 atoms for simplicity
        _write_itp(tmp_path, prefix, "water", WATER_ITP_CONTENT)

        # Write .top with molecule layout
        top_path = tmp_path / f"{prefix}.top"
        top_path.write_text("""\
[ system ]
test

[ molecules ]
MOL0    1
MOL1    2
water   3
""")

        # Build the broken .gro
        # Protein: 6 atoms, 2 residues (resids 1,1,1,2,2,2) — correct
        protein_lines = [
            _make_gro_line(1, "ALA", "N", 1),
            _make_gro_line(1, "ALA", "CA", 2),
            _make_gro_line(1, "ALA", "H", 3),
            _make_gro_line(2, "GLY", "N", 4),
            _make_gro_line(2, "GLY", "CA", 5),
            _make_gro_line(2, "GLY", "H", 6),
        ]

        # Polymer copy 0: 10 atoms, resids 1,1,1,2,2,2,2,3,3,3 — correct
        poly0_lines = [
            _make_gro_line(1, "RES1", "C1", 7),
            _make_gro_line(1, "RES1", "H1", 8),
            _make_gro_line(1, "RES1", "O1", 9),
            _make_gro_line(2, "RES2", "C2", 10),
            _make_gro_line(2, "RES2", "H2", 11),
            _make_gro_line(2, "RES2", "N1", 12),
            _make_gro_line(2, "RES2", "H3", 13),
            _make_gro_line(3, "RES3", "C3", 14),
            _make_gro_line(3, "RES3", "H4", 15),
            _make_gro_line(3, "RES3", "S1", 16),
        ]

        # Polymer copy 1: broken — resids shifted by +1 (copy_index=1)
        poly1_lines = [
            _make_gro_line(2, "RES1", "C1", 17),  # Should be 4
            _make_gro_line(2, "RES1", "H1", 18),
            _make_gro_line(2, "RES1", "O1", 19),
            _make_gro_line(3, "RES2", "C2", 20),  # Should be 5
            _make_gro_line(3, "RES2", "H2", 21),
            _make_gro_line(3, "RES2", "N1", 22),
            _make_gro_line(3, "RES2", "H3", 23),
            _make_gro_line(4, "RES3", "C3", 24),  # Should be 6
            _make_gro_line(4, "RES3", "H4", 25),
            _make_gro_line(4, "RES3", "S1", 26),
        ]

        # Water: 3 copies, single residue each — correct numbering
        water_lines = [
            _make_gro_line(1, "HOH", "OW", 27),
            _make_gro_line(1, "HOH", "HW1", 28),
            _make_gro_line(1, "HOH", "HW2", 29),
            _make_gro_line(2, "HOH", "OW", 30),
            _make_gro_line(2, "HOH", "HW1", 31),
            _make_gro_line(2, "HOH", "HW2", 32),
            _make_gro_line(3, "HOH", "OW", 33),
            _make_gro_line(3, "HOH", "HW1", 34),
            _make_gro_line(3, "HOH", "HW2", 35),
        ]

        all_atom_lines = protein_lines + poly0_lines + poly1_lines + water_lines
        n_atoms = len(all_atom_lines)
        gro_lines = (
            ["Generated by test"]
            + [f" {n_atoms}"]
            + all_atom_lines
            + ["   5.0000000   5.0000000   5.0000000"]
        )
        gro_path = tmp_path / f"{prefix}.gro"
        gro_path.write_text("\n".join(gro_lines) + "\n")

        return gro_path, top_path, prefix

    def test_fixes_polymer_residue_numbering(self, tmp_path):
        gro_path, top_path, prefix = self._setup_system(tmp_path)

        # Create a minimal GromacsExporter (no Interchange needed for this method)
        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro_path, top_path, tmp_path, prefix)

        # Parse fixed .gro
        lines = gro_path.read_text().splitlines()
        atom_lines = lines[2:-1]

        # Extract residue numbers
        resids = [int(line[:5].strip()) for line in atom_lines]

        # Protein (6 atoms, 2 residues): globally sequential 1, 1, 1, 2, 2, 2
        assert resids[0:6] == [1, 1, 1, 2, 2, 2], f"Protein resids: {resids[0:6]}"

        # Polymer copy 0 (10 atoms, 3 residues): continues from protein → 3, 3, 3, 4, 4, 4, 4, 5, 5, 5
        assert resids[6:16] == [
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
        ], f"Polymer copy 0 resids: {resids[6:16]}"

        # Polymer copy 1 (10 atoms, 3 residues): continues → 6, 6, 6, 7, 7, 7, 7, 8, 8, 8
        assert resids[16:26] == [
            6,
            6,
            6,
            7,
            7,
            7,
            7,
            8,
            8,
            8,
        ], f"Polymer copy 1 resids: {resids[16:26]}"

        # Water (9 atoms): single-residue, numbering preserved as-is
        assert resids[26:35] == [1, 1, 1, 2, 2, 2, 3, 3, 3], f"Water resids: {resids[26:35]}"

    def test_preserves_coordinates(self, tmp_path):
        gro_path, top_path, prefix = self._setup_system(tmp_path)

        # Read coordinates before fix
        lines_before = gro_path.read_text().splitlines()
        coords_before = [line[20:] for line in lines_before[2:-1]]

        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro_path, top_path, tmp_path, prefix)

        # Read coordinates after fix
        lines_after = gro_path.read_text().splitlines()
        coords_after = [line[20:] for line in lines_after[2:-1]]

        assert coords_before == coords_after

    def test_preserves_box_vectors(self, tmp_path):
        gro_path, top_path, prefix = self._setup_system(tmp_path)

        lines_before = gro_path.read_text().splitlines()
        box_before = lines_before[-1]

        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro_path, top_path, tmp_path, prefix)

        lines_after = gro_path.read_text().splitlines()
        box_after = lines_after[-1]

        assert box_before == box_after

    def test_preserves_atom_names(self, tmp_path):
        gro_path, top_path, prefix = self._setup_system(tmp_path)

        lines_before = gro_path.read_text().splitlines()
        names_before = [line[10:15] for line in lines_before[2:-1]]

        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro_path, top_path, tmp_path, prefix)

        lines_after = gro_path.read_text().splitlines()
        names_after = [line[10:15] for line in lines_after[2:-1]]

        assert names_before == names_after

    def test_no_op_for_single_residue_only(self, tmp_path):
        """If all molecules are single-residue, .gro should be unchanged."""
        prefix = "sys"
        _write_itp(tmp_path, prefix, "water", WATER_ITP_CONTENT)

        top = tmp_path / f"{prefix}.top"
        top.write_text("[ system ]\ntest\n[ molecules ]\nwater   2\n")

        gro = tmp_path / f"{prefix}.gro"
        gro_content = (
            "\n".join(
                [
                    "title",
                    "6",
                    _make_gro_line(1, "HOH", "OW", 1),
                    _make_gro_line(1, "HOH", "HW1", 2),
                    _make_gro_line(1, "HOH", "HW2", 3),
                    _make_gro_line(2, "HOH", "OW", 4),
                    _make_gro_line(2, "HOH", "HW1", 5),
                    _make_gro_line(2, "HOH", "HW2", 6),
                    "   3.0000000   3.0000000   3.0000000",
                ]
            )
            + "\n"
        )
        gro.write_text(gro_content)

        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro, top, tmp_path, prefix)

        # Content should be unchanged (no multi-residue molecules)
        assert gro.read_text() == gro_content


class TestFixGroMultiplePolymerTypes:
    """Test residue numbering with multiple distinct polymer molecule types."""

    def test_two_polymer_types(self, tmp_path):
        """Two different polymer types, each with 2 copies.

        Layout: protein(1) + polyA(2 copies) + polyB(2 copies)
        polyA: 10 atoms, 3 residues
        polyB: 10 atoms, 3 residues (different sequence, same atom count)
        """
        prefix = "sys"

        # Protein
        _write_itp(tmp_path, prefix, "MOL0", PROTEIN_ITP_CONTENT)

        # Two polymer types
        _write_itp(
            tmp_path,
            prefix,
            "MOL1",
            POLYMER_ITP_TEMPLATE.format(mol_name="MOL1"),
        )
        _write_itp(
            tmp_path,
            prefix,
            "MOL2",
            POLYMER_ITP_TEMPLATE.format(mol_name="MOL2"),
        )

        top = tmp_path / f"{prefix}.top"
        top.write_text("[ system ]\ntest\n[ molecules ]\nMOL0    1\nMOL1    2\nMOL2    2\n")

        # Build GRO with broken numbering
        lines = ["title", "46"]

        # Protein: 6 atoms
        for i, (resid, rname, aname) in enumerate(
            [
                (1, "ALA", "N"),
                (1, "ALA", "CA"),
                (1, "ALA", "H"),
                (2, "GLY", "N"),
                (2, "GLY", "CA"),
                (2, "GLY", "H"),
            ],
            start=1,
        ):
            lines.append(_make_gro_line(resid, rname, aname, i))

        atom_num = 7
        # polyA copy 0: resids 1,1,1,2,2,2,2,3,3,3
        for resid, rname, aname in [
            (1, "RES1", "C1"),
            (1, "RES1", "H1"),
            (1, "RES1", "O1"),
            (2, "RES2", "C2"),
            (2, "RES2", "H2"),
            (2, "RES2", "N1"),
            (2, "RES2", "H3"),
            (3, "RES3", "C3"),
            (3, "RES3", "H4"),
            (3, "RES3", "S1"),
        ]:
            lines.append(_make_gro_line(resid, rname, aname, atom_num))
            atom_num += 1

        # polyA copy 1: broken (resids +1)
        for resid, rname, aname in [
            (2, "RES1", "C1"),
            (2, "RES1", "H1"),
            (2, "RES1", "O1"),
            (3, "RES2", "C2"),
            (3, "RES2", "H2"),
            (3, "RES2", "N1"),
            (3, "RES2", "H3"),
            (4, "RES3", "C3"),
            (4, "RES3", "H4"),
            (4, "RES3", "S1"),
        ]:
            lines.append(_make_gro_line(resid, rname, aname, atom_num))
            atom_num += 1

        # polyB copy 0: resids 1,1,1,2,2,2,2,3,3,3
        for resid, rname, aname in [
            (1, "RES1", "C1"),
            (1, "RES1", "H1"),
            (1, "RES1", "O1"),
            (2, "RES2", "C2"),
            (2, "RES2", "H2"),
            (2, "RES2", "N1"),
            (2, "RES2", "H3"),
            (3, "RES3", "C3"),
            (3, "RES3", "H4"),
            (3, "RES3", "S1"),
        ]:
            lines.append(_make_gro_line(resid, rname, aname, atom_num))
            atom_num += 1

        # polyB copy 1: broken (resids +1)
        for resid, rname, aname in [
            (2, "RES1", "C1"),
            (2, "RES1", "H1"),
            (2, "RES1", "O1"),
            (3, "RES2", "C2"),
            (3, "RES2", "H2"),
            (3, "RES2", "N1"),
            (3, "RES2", "H3"),
            (4, "RES3", "C3"),
            (4, "RES3", "H4"),
            (4, "RES3", "S1"),
        ]:
            lines.append(_make_gro_line(resid, rname, aname, atom_num))
            atom_num += 1

        lines.append("   5.0000000   5.0000000   5.0000000")
        gro_path = tmp_path / f"{prefix}.gro"
        gro_path.write_text("\n".join(lines) + "\n")

        # Fix
        exporter = GromacsExporter.__new__(GromacsExporter)
        exporter._fix_gro_residue_numbering(gro_path, top, tmp_path, prefix)

        # Parse fixed
        fixed_lines = gro_path.read_text().splitlines()
        resids = [int(line[:5].strip()) for line in fixed_lines[2:-1]]

        # Protein: 1,1,1,2,2,2 (globally sequential)
        assert resids[0:6] == [1, 1, 1, 2, 2, 2]

        # polyA copy 0: 3,3,3,4,4,4,4,5,5,5 (continues from protein resid 2)
        assert resids[6:16] == [3, 3, 3, 4, 4, 4, 4, 5, 5, 5]

        # polyA copy 1: 6,6,6,7,7,7,7,8,8,8 (continues from polyA copy 0 resid 5)
        assert resids[16:26] == [6, 6, 6, 7, 7, 7, 7, 8, 8, 8]

        # polyB copy 0: 9,9,9,10,10,10,10,11,11,11 (continues from polyA copy 1 resid 8)
        assert resids[26:36] == [9, 9, 9, 10, 10, 10, 10, 11, 11, 11]

        # polyB copy 1: 12,12,12,13,13,13,13,14,14,14
        assert resids[36:46] == [12, 12, 12, 13, 13, 13, 13, 14, 14, 14]


class TestEnergyMinimizationHelpers:
    """Tests for EM health checks."""

    def test_run_script_template_resource_exists(self):
        """Local GROMACS run template should be packaged as a resource."""
        template = resources.files("polyzymd.engines.gromacs").joinpath(
            "templates", "run_gromacs.sh.jinja"
        )

        assert template.is_file()
        assert "GROMACS Workflow Script" in template.read_text()

    def test_run_script_generator_uses_single_stage_em_with_health_check(self, tmp_path):
        """Generated local run script should use single-stage EM and preserve health check."""
        script_path = tmp_path / "run_test_gromacs.sh"
        generator = RunScriptGenerator(prefix="system", equilibration_mdps=["eq_01_nvt.mdp"])

        generator.generate(script_path)
        script_content = script_path.read_text()

        assert "em.mdp" in script_content
        assert "em.gro" in script_content
        assert "em.tpr" in script_content

        assert "em_soft.mdp" not in script_content
        assert "em_soft.gro" not in script_content
        assert "em_soft.tpr" not in script_content

        assert "grep -qi" in script_content
        assert "force.*not finite" in script_content
        _assert_no_unresolved_jinja(script_content)

    def test_run_script_generator_preserves_command_semantics(self, tmp_path):
        """Generated local run script should preserve GROMACS command sequence."""
        script_path = tmp_path / "run_test_gromacs.sh"
        generator = RunScriptGenerator(
            prefix="system",
            equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
            gmx_command="/opt/gromacs/bin/gmx",
        )

        generator.generate(script_path)
        script_content = script_path.read_text()

        assert 'GMX="/opt/gromacs/bin/gmx"' in script_content
        assert 'PREFIX="system"' in script_content
        assert (
            "$GMX grompp -f em.mdp -c ${PREFIX}.gro -r ${PREFIX}.gro "
            "-p ${PREFIX}.top -o em.tpr -maxwarn 1"
        ) in script_content
        assert "$GMX mdrun -deffnm em -v" in script_content
        assert "$GMX grompp -f eq_01_nvt.mdp -c em.gro -r em.gro" in script_content
        assert "$GMX mdrun -deffnm eq_01 -v" in script_content
        assert "$GMX grompp -f eq_02_npt.mdp -c eq_01.gro -r em.gro -t eq_01.cpt" in script_content
        assert "$GMX mdrun -deffnm eq_02 -v" in script_content
        assert "$GMX grompp -f prod.mdp -c ${LAST_EQ}.gro -t ${LAST_EQ}.cpt" in script_content
        assert "$GMX mdrun -deffnm prod -v" in script_content
        assert "prod_nojump.xtc" in script_content
        assert "prod_centered.xtc" in script_content

    @pytest.mark.parametrize(
        ("field", "kwargs"),
        [
            ("gmx_command", {"gmx_command": "/opt/gromacs/bin/gmx$1"}),
            ("gmx_command", {"gmx_command": r"C:\gromacs\bin\gmx"}),
            ("gmx_command", {"gmx_command": "gmx mpi"}),
            ("prefix", {"prefix": "bad$prefix"}),
            ("prefix", {"prefix": r"bad\prefix"}),
            ("prefix", {"prefix": "bad prefix"}),
            ("equilibration_mdp", {"equilibration_mdps": ["eq_$1.mdp"]}),
            ("equilibration_mdp", {"equilibration_mdps": [r"eq\01.mdp"]}),
            ("equilibration_mdp", {"equilibration_mdps": ["eq 01.mdp"]}),
        ],
    )
    def test_run_script_generator_rejects_unsafe_shell_tokens(self, tmp_path, field, kwargs):
        """Local run script values should reject unsafe shell tokens."""
        script_path = tmp_path / "run_test_gromacs.sh"
        init_kwargs = {
            "prefix": "system",
            "equilibration_mdps": ["eq_01_nvt.mdp"],
            "gmx_command": "gmx",
        }
        init_kwargs.update(kwargs)
        generator = RunScriptGenerator(**init_kwargs)

        with pytest.raises(ValueError, match=field):
            generator.generate(script_path)

    def test_em_health_check_in_runner(self, tmp_path):
        """Runner health check should fail on non-finite force signatures."""
        em_log = tmp_path / "em.log"
        em_log.write_text("force on at least one atom is not finite\n")

        runner = GromacsRunner(
            working_dir=tmp_path,
            prefix="system",
            equilibration_mdps=[],
        )

        with pytest.raises(RuntimeError, match="infinite forces"):
            runner._check_energy_minimization_health()
