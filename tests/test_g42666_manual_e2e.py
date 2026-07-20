"""Tests for the guarded G42666 manual E2E support script."""

from __future__ import annotations

import math

import pytest

from tests._support.g42666_manual_e2e import _read_state_rows


def test_read_state_rows_preserves_commented_openmm_header(tmp_path):
    """Parse OpenMM StateDataReporter CSV output with its commented header."""
    state_path = tmp_path / "state.csv"
    rows = [
        f'{step},"{energy:.6f}","{temperature:.3f}","{density:.6f}"'
        for step, energy, temperature, density in (
            (index * 1000, -1250.0 + index * 0.25, 300.0 + index * 0.01, 0.998 + index * 1e-5)
            for index in range(1, 151)
        )
    ]
    state_path.write_text(
        "# OpenMM StateDataReporter metadata\n"
        '#"Step","Potential Energy (kJ/mole)","Temperature (K)","Density (g/mL)"\n'
        + "\n".join(rows)
        + "\n# trailing diagnostic comment\n",
        encoding="utf-8",
    )

    parsed = _read_state_rows(state_path)

    assert len(parsed) == 150
    assert math.isclose(parsed[0]["potential_energy"], -1249.75)
    assert math.isclose(parsed[-1]["potential_energy"], -1212.5)


def test_read_state_rows_ignores_misleading_comment_with_header_words(tmp_path):
    """Ignore non-CSV comments that mention StateDataReporter column names."""
    state_path = tmp_path / "state.csv"
    state_path.write_text(
        "# Step Potential Energy is described in this metadata note\n"
        '#"Step","Potential Energy (kJ/mole)","Temperature (K)"\n'
        '1000,"-1.500000","300.000"\n',
        encoding="utf-8",
    )

    parsed = _read_state_rows(state_path)

    assert parsed == [{"potential_energy": -1.5}]


def test_read_state_rows_rejects_duplicate_real_headers(tmp_path):
    """Reject a second structurally valid StateDataReporter header."""
    state_path = tmp_path / "state.csv"
    state_path.write_text(
        '#"Step","Potential Energy (kJ/mole)","Temperature (K)"\n'
        '1000,"-1.500000","300.000"\n'
        '#"Step","Potential Energy (kJ/mole)","Temperature (K)"\n'
        '2000,"-1.250000","300.100"\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="multiple StateDataReporter headers"):
        _read_state_rows(state_path)
