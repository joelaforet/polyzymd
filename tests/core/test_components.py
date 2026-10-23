"""Component records: which atoms a component owns and where their parameters come from."""

from __future__ import annotations

import pytest

from polyzymd.core.components import (
    DEFAULT_SEGIDS,
    AtomRange,
    ChargeSource,
    ComponentRecord,
    Domain,
    check_layout,
    merge_ranges,
    molecule_ranges,
    validate_segid,
)


def _record(name="protein", segid="PROT", start=0, stop=10, domains=None):
    atoms = [AtomRange(start, stop)]
    domains = domains or [Domain(name, tuple(atoms), ("ff.offxml",), ChargeSource("library"))]
    return ComponentRecord(name, "protein", segid, "fixed", atoms, domains)


@pytest.mark.parametrize("segid", ["", "TOOLONG", "A B", "ÉTÉ"])
def test_segment_ids_must_fit_the_psf(segid) -> None:
    """PSF segment IDs hold at most four plain characters."""
    with pytest.raises(ValueError, match="segment ID"):
        validate_segid(segid)


def test_a_reserved_segment_id_is_refused() -> None:
    """A user-defined component cannot take a built-in type's segment ID."""
    with pytest.raises(ValueError, match="reserved"):
        validate_segid("PROT", reserved=DEFAULT_SEGIDS.values())


def test_domains_must_cover_each_atom_once() -> None:
    """A component whose domains miss or repeat atoms is refused."""
    half = Domain("a", (AtomRange(0, 5),), ("ff",), ChargeSource("library"))
    with pytest.raises(ValueError, match="exactly once"):
        _record(domains=[half])
    overlap = Domain("b", (AtomRange(3, 10),), ("ff",), ChargeSource("library"))
    with pytest.raises(ValueError, match="exactly once"):
        _record(domains=[half, overlap])


def test_a_modified_protein_is_two_domains() -> None:
    """The 1.4 shape: one component, protein and site domains with different sources."""
    protein = Domain("protein", (AtomRange(0, 7),), ("ff14sb.offxml",), ChargeSource("library"))
    site = Domain(
        "site", (AtomRange(7, 10),), ("openff-2.3.0.offxml",), ChargeSource("nagl", "model.pt")
    )
    record = _record(domains=[protein, site])
    data = record.to_dict()
    assert [d["name"] for d in data["domains"]] == ["protein", "site"]
    assert data["domains"][1]["charges"] == {
        "method": "nagl",
        "model": "model.pt",
        "source_file": None,
    }
    assert "later domain" in data["boundary_rule"]


def test_layout_must_cover_the_system_once_with_unique_handles() -> None:
    """Components must tile the system's atoms, with distinct names and segment IDs."""
    first = _record()
    second = _record(name="ligand", segid="LIG", start=10, stop=12)
    check_layout([first, second], 12)
    with pytest.raises(ValueError, match="cover 12 of 13"):
        check_layout([first, second], 13)
    with pytest.raises(ValueError, match="segid 'PROT'"):
        check_layout([first, _record(name="other", start=10, stop=12)], 12)


def test_adjacent_molecules_merge_into_one_range() -> None:
    """Consecutive molecules are written as one range in the record."""
    ranges = molecule_ranges([3, 3, 4], start=10)
    assert ranges == [AtomRange(10, 13), AtomRange(13, 16), AtomRange(16, 20)]
    assert merge_ranges(ranges) == [AtomRange(10, 20)]
