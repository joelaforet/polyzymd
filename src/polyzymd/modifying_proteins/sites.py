"""Residue sites written as ``A:N58``, ``A:ASN58`` or ``A:N58A``.

A site names a chain, a residue code and a residue number, with an optional
insertion code. The code is checked against the structure, so a numbering
offset fails with a message instead of modifying the wrong residue.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: One-letter codes of the canonical amino acids.
ONE_LETTER = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

#: The canonical residues: mutating to one of these keeps ff14SB parameters.
CANONICAL = frozenset(ONE_LETTER.values())

#: Protonation and disulfide variants that name the same canonical residue.
_VARIANTS = {
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYM": "CYS",
    "CYX": "CYS",
}

_SITE = re.compile(r"^(?P<chain>[A-Za-z0-9]+):(?P<code>.+?)(?P<number>-?\d+)(?P<icode>[A-Za-z]?)$")
_RESIDUE = re.compile(r"^(?P<code>.+?)(?P<number>-?\d+)(?P<icode>[A-Za-z]?)$")


@dataclass(frozen=True)
class Site:
    """One residue of a protein: chain, residue code, number, insertion code."""

    chain_id: str
    code: str
    resnum: int
    icode: str = ""

    def __str__(self) -> str:
        return f"{self.chain_id}:{self.code}{self.resnum}{self.icode}"

    @property
    def resname(self) -> str:
        """The three-letter (or CCD) residue name the code stands for."""
        return ONE_LETTER.get(self.code, self.code)


def parse_site(text: str) -> Site:
    """Parse ``A:N58``, ``A:ASN58``, ``A:N58A`` or ``AB:4II12`` into a Site."""
    match = _SITE.match(str(text).strip())
    if not match or not _valid_code(match["code"]):
        raise ValueError(
            f"Cannot read site {text!r}. Write chain:residue, for example 'A:N58', "
            "'A:ASN58' or, with an insertion code, 'A:N58A'."
        )
    code = match["code"].upper()
    return Site(match["chain"], code, int(match["number"]), match["icode"].upper())


def expand_sites(chains: list[str], residues: list[str]) -> list[Site]:
    """Return every residue of ``residues`` (``N58``) on every chain of ``chains``."""
    sites = []
    for chain in chains:
        for residue in residues:
            match = _RESIDUE.match(str(residue).strip())
            if not match or not _valid_code(match["code"]):
                raise ValueError(
                    f"Cannot read residue {residue!r}; with `chains`, write residues "
                    "without a chain, for example 'N58' or 'ASN58'."
                )
            sites.append(
                Site(chain, match["code"].upper(), int(match["number"]), match["icode"].upper())
            )
    return sites


def _valid_code(code: str) -> bool:
    if not code.isalnum():
        return False
    return code.upper() in ONE_LETTER if len(code) == 1 else True


def canonical_name(resname: str) -> str:
    """Return the canonical residue a variant name stands for (``LYN`` -> ``LYS``)."""
    return _VARIANTS.get(resname.upper(), resname.upper())


def resolve(protein, site: Site):
    """Return the mBuild residue at ``site``, checking its name.

    Raises ValueError when the residue does not exist or has another name,
    for example after a numbering offset.
    """
    from mbuild.exceptions import MBuildError

    try:
        residue = protein.get_residue(site.resnum, chain_id=site.chain_id, icode=site.icode)
    except MBuildError as error:
        raise ValueError(f"Site {site}: {error}") from error
    found = residue.name.upper()
    if canonical_name(found) != canonical_name(site.resname):
        raise ValueError(
            f"Site {site} names {site.resname}, but residue {site.resnum}{site.icode} "
            f"of chain {site.chain_id} is {found}. Check the residue numbering of the "
            "structure, or write the residue it holds now (after a mutation, its new name)."
        )
    return residue
