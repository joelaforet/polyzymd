"""Reading residue sites and checking them against a structure."""

from pathlib import Path

import pytest

from polyzymd.modifying_proteins.sites import Site, expand_sites, parse_site, resolve

TRPCAGE = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("A:N58", Site("A", "N", 58)),
        ("A:ASN58", Site("A", "ASN", 58)),
        ("A:N58A", Site("A", "N", 58, "A")),
        ("AB:4II12", Site("AB", "4II", 12)),
        ("a:k8", Site("a", "K", 8)),
    ],
)
def test_parse_site_reads_one_and_three_letter_codes(text, expected) -> None:
    assert parse_site(text) == expected


@pytest.mark.parametrize("text", ["N58", "A:58", "A:X58", "A:N", "A:N-5x1"])
def test_parse_site_rejects_unreadable_sites(text) -> None:
    with pytest.raises(ValueError, match="Cannot read site"):
        parse_site(text)


def test_site_names_its_residue() -> None:
    assert parse_site("A:N58").resname == "ASN"
    assert parse_site("A:4II12").resname == "4II"
    assert str(parse_site("A:ASN58A")) == "A:ASN58A"


def test_expand_sites_puts_each_residue_on_each_chain() -> None:
    sites = expand_sites(["A", "C"], ["N58", "ASN107"])
    assert [str(site) for site in sites] == ["A:N58", "A:ASN107", "C:N58", "C:ASN107"]


def test_expand_sites_rejects_a_residue_with_a_chain() -> None:
    with pytest.raises(ValueError, match="without a chain"):
        expand_sites(["A"], ["A:N58"])


def test_resolve_checks_the_residue_name() -> None:
    from mbuild.biopolymers import Protein

    protein = Protein(str(TRPCAGE))
    assert resolve(protein, parse_site("A:K8")).name == "LYS"
    assert resolve(protein, parse_site("A:LYS8")).name == "LYS"
    # One residue off: residue 9 is Asp, so a numbering offset fails loudly.
    with pytest.raises(ValueError, match="residue 9 of chain A is ASP"):
        resolve(protein, parse_site("A:K9"))
    with pytest.raises(ValueError, match="No residue with number 99"):
        resolve(protein, parse_site("A:K99"))


class _Residues:
    """A stand-in structure with residues named by (number, name)."""

    def __init__(self, residues):
        self.residues = residues

    def get_residue(self, resnum, chain_id=None, icode=""):
        from mbuild.exceptions import MBuildError

        if resnum not in self.residues:
            raise MBuildError(f"No residue with number {resnum}")
        return type("Residue", (), {"name": self.residues[resnum], "resnum": resnum})()


def test_resolve_text_reads_a_code_ending_in_a_digit_by_the_structure() -> None:
    from polyzymd.modifying_proteins.sites import candidate_sites, resolve_text

    assert [(s.code, s.resnum) for s in candidate_sites("A:HY312")] == [
        ("HY", 312),
        ("HY3", 12),
        ("HY31", 2),
    ]
    site, residue = resolve_text(_Residues({12: "HY3"}), "A:HY312")
    assert (site.code, site.resnum, residue.name) == ("HY3", 12, "HY3")


def test_resolve_text_refuses_two_readings_that_both_exist() -> None:
    from polyzymd.modifying_proteins.sites import resolve_text

    with pytest.raises(ValueError, match="more than one residue"):
        resolve_text(_Residues({12: "HY3", 312: "HY"}), "A:HY312")


def test_resolve_text_keeps_the_numbering_hint_when_nothing_matches() -> None:
    from polyzymd.modifying_proteins.sites import resolve_text

    with pytest.raises(ValueError, match="residue 9 of chain A is ASP"):
        resolve_text(_Residues({9: "ASP"}), "A:K9")
    with pytest.raises(ValueError, match="Other readings tried"):
        resolve_text(_Residues({312: "GLY"}), "A:HY312")
