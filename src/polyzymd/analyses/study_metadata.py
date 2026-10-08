"""The publishing metadata of a study: ``metadata:`` in ``study.yaml``.

:func:`check_metadata` checks the block and fills what is missing with TODO
placeholders, returning a warning for each, because publishing metadata
must never stop an analysis or a freeze. :func:`citation_cff` and
:func:`zenodo_json` write ``CITATION.cff`` (Citation File Format 1.2.0) and
``.zenodo.json`` from the one block, so the two cannot disagree; Zenodo
reads only ``.zenodo.json`` when both exist. Both cite PolyzyMD
(:mod:`polyzymd.citation`). The fields follow the FAIR principles
(Wilkinson et al. 2016), and the minimum metadata of Amaro et al. (2025),
which asks in particular for the purpose of the simulations; see the
"Study folders" explanation page.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

TODO = "TODO"
_KEYS = (
    "doi",
    "title",
    "description",
    "purpose",
    "keywords",
    "system_type",
    "authors",
    "contact",
    "license",
    "funding",
    "related",
    "zenodo",
)
_PERSON = ("family-names", "given-names", "name-suffix", "orcid", "affiliation", "email")
_LICENSE = ("data", "code")
_RELATED = ("paper", "trajectories", "experimental")
_PAPER = ("title", "doi", "status", "journal", "year")
_TRAJECTORIES = ("doi", "conditions", "title")
_EXPERIMENTAL = ("doi", "description")
_FUNDING = ("funder", "award", "funder_doi")
_ZENODO = ("communities", "access_right")
#: Default licences, matching the LICENSE files polyzymd study init writes.
DEFAULT_LICENSE = {"data": "CC-BY-4.0", "code": "MIT"}
#: Values of ``related.paper.status``, those of Citation File Format 1.2.0.
PAPER_STATUS = ("in-preparation", "submitted", "in-press", "preprint", "advance-online", "abstract")
#: Values of ``zenodo.access_right``.
ACCESS = ("open", "embargoed", "restricted", "closed")


#: Placeholder parts of a DOI: ``NNNN`` or ``XXXX`` runs, ``...``, or a last part of zeros.
_PLACEHOLDER_DOI = re.compile(r"N{4,}|X{4,}|\.\.\.|[./]0+$", re.IGNORECASE)


def is_placeholder(doi: Any) -> bool:
    """Return whether ``doi`` is missing or still a placeholder such as ``10.5281/zenodo.NNNNNNN``.

    Placeholders are ``XXXX`` or ``NNNN`` runs, ``...``, a last part of only
    zeros (``zenodo.0000000``), ``placeholder`` and ``TODO``.
    """
    text = str(doi or "").strip()
    return (
        not text
        or bool(_PLACEHOLDER_DOI.search(text))
        or "placeholder" in text.lower()
        or TODO in text
    )


def _known(raw: Mapping, keys: tuple[str, ...], where: str) -> None:
    from polyzymd.analyses.study_file import _unknown

    _unknown(raw, keys, where)


def _mapping(value: Any, where: str) -> Mapping:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ProtocolError(f"{where} must be a mapping, got {value!r}.", hint="See study_folders.")
    return value


def _sequence(value: Any, where: str) -> list:
    if value is None:
        return []
    if isinstance(value, (str, Mapping)) or not isinstance(value, (list, tuple)):
        raise ProtocolError(f"{where} must be a list, got {value!r}.", hint="Write it as [a, b].")
    return list(value)


def check_metadata(raw: Any, what: str = "study") -> tuple[dict[str, Any], list[str]]:
    """Check the ``metadata:`` block and return it completed, with a warning per gap.

    Missing ``title``, ``description``, ``purpose`` and ``authors`` become
    ``TODO`` placeholders; a missing licence becomes the default
    (:data:`DEFAULT_LICENSE`); a missing or placeholder paper DOI is noted.
    ``doi`` is the dataset's own DOI, reserved in Zenodo before publishing;
    it is ``None`` until set, with a warning that names ``what`` is published
    (``"study"``, or ``"project"`` for a project's metadata).

    Raises
    ------
    ProtocolError
        For an unknown key (with the nearest known one) or a value of the
        wrong type. Gaps never raise.
    """
    raw = _mapping(raw, "metadata")
    _known(raw, _KEYS, "metadata")
    warnings: list[str] = []
    meta: dict[str, Any] = {}
    for key in ("title", "description", "purpose"):
        value = raw.get(key)
        if not value or not str(value).strip():
            warnings.append(f"metadata.{key} is missing")
            meta[key] = f"{TODO}: add metadata.{key} to {what}.yaml"
        else:
            meta[key] = str(value).strip()
    for key in ("keywords", "system_type"):
        meta[key] = [str(item) for item in _sequence(raw.get(key), f"metadata.{key}")]
    if not meta["keywords"]:
        warnings.append("metadata.keywords is empty; keywords make the deposit findable")

    authors = []
    for index, person in enumerate(_sequence(raw.get("authors"), "metadata.authors")):
        person = _mapping(person, f"metadata.authors[{index}]")
        _known(person, _PERSON, f"metadata.authors[{index}]")
        if not person.get("family-names"):
            raise ProtocolError(
                f"metadata.authors[{index}] has no family-names.",
                hint="Write each author as {family-names: ..., given-names: ..., orcid: ...}.",
            )
        if not person.get("orcid"):
            warnings.append(f"author {person['family-names']} has no orcid")
        authors.append({k: str(v) for k, v in person.items() if v})
    if not authors:
        warnings.append("metadata.authors is empty")
        authors = [{"name": f"{TODO}: add metadata.authors to {what}.yaml"}]
    meta["authors"] = authors
    meta["contact"] = dict(_mapping(raw.get("contact"), "metadata.contact"))
    doi = raw.get("doi")
    if is_placeholder(doi):
        warnings.append(
            f"metadata.doi is not set: reserve a DOI for the {what} in Zenodo, add it here and "
            "refreeze"
        )
        meta["doi"] = None
    else:
        meta["doi"] = str(doi).strip()

    license_ = _mapping(raw.get("license"), "metadata.license")
    _known(license_, _LICENSE, "metadata.license")
    meta["license"] = {**DEFAULT_LICENSE, **{k: str(v) for k, v in license_.items()}}

    funding = []
    for index, grant in enumerate(_sequence(raw.get("funding"), "metadata.funding")):
        grant = _mapping(grant, f"metadata.funding[{index}]")
        _known(grant, _FUNDING, f"metadata.funding[{index}]")
        grant = {k: str(v) for k, v in grant.items() if v}
        if "funder_doi" in grant and is_placeholder(grant["funder_doi"]):
            warnings.append(f"metadata.funding[{index}].funder_doi is a placeholder; left out")
            del grant["funder_doi"]
        funding.append(grant)
    meta["funding"] = funding

    related = _mapping(raw.get("related"), "metadata.related")
    _known(related, _RELATED, "metadata.related")
    paper = dict(_mapping(related.get("paper"), "metadata.related.paper"))
    _known(paper, _PAPER, "metadata.related.paper")
    if paper.get("status") and paper["status"] not in PAPER_STATUS:
        raise ProtocolError(
            f"metadata.related.paper.status {paper['status']!r} is not one of {', '.join(PAPER_STATUS)}.",
            hint="Leave status out once the paper is published.",
        )
    if is_placeholder(paper.get("doi")):
        warnings.append("the paper DOI is missing or a placeholder; refreeze once it is known")
    trajectories = []
    for index, entry in enumerate(
        _sequence(related.get("trajectories"), "metadata.related.trajectories")
    ):
        entry = _mapping(entry, f"metadata.related.trajectories[{index}]")
        _known(entry, _TRAJECTORIES, f"metadata.related.trajectories[{index}]")
        if is_placeholder(entry.get("doi")):
            warnings.append(f"trajectory deposit {index + 1} has no DOI yet")
        trajectories.append(
            {
                **{k: v for k, v in entry.items() if k != "conditions"},
                "conditions": [str(c) for c in _sequence(entry.get("conditions"), "conditions")],
            }
        )
    if not trajectories:
        warnings.append(
            "metadata.related.trajectories lists no trajectory deposit; reproducing the analyses "
            "needs the trajectories"
        )
    experimental = []
    for index, entry in enumerate(
        _sequence(related.get("experimental"), "metadata.related.experimental")
    ):
        entry = _mapping(entry, f"metadata.related.experimental[{index}]")
        _known(entry, _EXPERIMENTAL, f"metadata.related.experimental[{index}]")
        experimental.append(dict(entry))
    meta["related"] = {"paper": paper, "trajectories": trajectories, "experimental": experimental}

    zenodo = dict(_mapping(raw.get("zenodo"), "metadata.zenodo"))
    _known(zenodo, _ZENODO, "metadata.zenodo")
    access = zenodo.get("access_right", "open")
    if access not in ACCESS:
        raise ProtocolError(
            f"metadata.zenodo.access_right {access!r} is not one of {', '.join(ACCESS)}.",
            hint="Use open unless the data must be embargoed.",
        )
    meta["zenodo"] = {
        "communities": [str(c) for c in _sequence(zenodo.get("communities"), "communities")],
        "access_right": access,
    }
    todo = _todo_fields(raw, "metadata")
    if todo:
        warnings.append(f"TODO placeholders left in {', '.join(todo)}")
    return meta, warnings


def _todo_fields(value: Any, where: str) -> list[str]:
    """Return the fields under ``where`` still holding a TODO placeholder, DOIs left out."""
    if isinstance(value, Mapping):
        return [f for key, item in value.items() for f in _todo_fields(item, f"{where}.{key}")]
    if isinstance(value, (list, tuple)):
        return [f for i, item in enumerate(value) for f in _todo_fields(item, f"{where}[{i}]")]
    # A placeholder DOI has its own warning.
    return [where] if TODO in str(value) and not where.endswith("doi") else []


def _orcid_id(orcid: str) -> str:
    """Return the bare ORCID iD of ``orcid``, which may be a URL."""
    match = re.search(r"(\d{4}-\d{4}-\d{4}-\d{3}[\dX])", orcid)
    return match.group(1) if match else orcid


def _cff_person(person: dict[str, str]) -> dict[str, str]:
    out = dict(person)
    if "orcid" in out and not out["orcid"].startswith("http"):
        out["orcid"] = f"https://orcid.org/{_orcid_id(out['orcid'])}"
    return out


def citation_cff(
    meta: dict[str, Any], *, version: str, released: str, commit: str | None
) -> dict[str, Any]:
    """Return the study's CITATION.cff, as a mapping to write as YAML.

    The study's paper is the ``preferred-citation``, because CFF asks that it
    be cited instead of the dataset; PolyzyMD and its paper go in
    ``references``, the works the study builds on, with the trajectory
    deposits.
    """
    from polyzymd.citation import paper_reference, software_reference

    cff: dict[str, Any] = {
        "cff-version": "1.2.0",
        "type": "dataset",
        "message": "If you use this study, please cite this dataset and PolyzyMD, "
        "which produced its analyses.",
        "title": meta["title"],
        "abstract": meta["description"],
        "authors": [_cff_person(p) for p in meta["authors"]],
        "version": version,
        "date-released": released,
        "license": meta["license"]["data"],
        "keywords": meta["keywords"] or ["molecular dynamics"],
    }
    if meta.get("doi"):
        cff["doi"] = meta["doi"]
    if commit:
        cff["commit"] = commit
    paper = meta["related"]["paper"]
    if paper.get("title") or paper.get("doi"):
        preferred: dict[str, Any] = {
            "type": "article",
            "authors": [dict(person) for person in cff["authors"]],
            "title": paper.get("title") or meta["title"],
        }
        for key in ("doi", "status", "journal", "year"):
            if paper.get(key) and not (key == "doi" and is_placeholder(paper[key])):
                preferred[key] = paper[key]
        cff["preferred-citation"] = preferred
        cff["message"] = (
            "If you use this study, please cite the article in preferred-citation, "
            "this dataset, and PolyzyMD, which produced its analyses."
        )
    references = [software_reference()]
    if paper_reference():
        references.append(paper_reference())
    for entry in meta["related"]["trajectories"]:
        reference: dict[str, Any] = {
            "type": "data",
            "authors": [dict(person) for person in cff["authors"]],
            "title": entry.get("title")
            or f"Trajectories of {', '.join(entry['conditions']) or 'the study'}",
        }
        if not is_placeholder(entry.get("doi")):
            reference["doi"] = entry["doi"]
        references.append(reference)
    cff["references"] = references
    return cff


def zenodo_json(
    meta: dict[str, Any], *, version: str, released: str, method: str
) -> dict[str, Any]:
    """Return the study's ``.zenodo.json``: Zenodo deposit metadata, citing PolyzyMD."""
    from polyzymd.citation import DOI, REPOSITORY, paper_reference

    creators = []
    for person in meta["authors"]:
        name = person.get("name") or ", ".join(
            part for part in (person.get("family-names"), person.get("given-names")) if part
        )
        creator: dict[str, str] = {"name": name}
        if person.get("affiliation"):
            creator["affiliation"] = person["affiliation"]
        if person.get("orcid"):
            creator["orcid"] = _orcid_id(person["orcid"])
        creators.append(creator)
    related = []
    paper = meta["related"]["paper"]
    if not is_placeholder(paper.get("doi")):
        related.append(
            {
                "identifier": paper["doi"],
                "relation": "isSupplementTo",
                "resource_type": "publication-article",
            }
        )
    related.append(
        {"identifier": DOI or REPOSITORY, "relation": "requires", "resource_type": "software"}
    )
    polyzymd_paper = paper_reference()
    if polyzymd_paper and polyzymd_paper.get("doi"):
        related.append(
            {
                "identifier": polyzymd_paper["doi"],
                "relation": "references",
                "resource_type": "publication-article",
            }
        )
    for entry in meta["related"]["trajectories"]:
        if not is_placeholder(entry.get("doi")):
            related.append(
                {"identifier": entry["doi"], "relation": "references", "resource_type": "dataset"}
            )
    for entry in meta["related"]["experimental"]:
        if entry.get("doi"):
            related.append({"identifier": entry["doi"], "relation": "references"})
    zenodo: dict[str, Any] = {
        "upload_type": "dataset",
        "title": meta["title"],
        "description": f"{meta['description']}\n\nPurpose: {meta['purpose']}",
        "creators": creators,
        "access_right": meta["zenodo"]["access_right"],
        "license": meta["license"]["data"].lower(),
        "version": version,
        "publication_date": released,
        **({"doi": meta["doi"]} if meta.get("doi") else {}),
        "keywords": meta["keywords"],
        "method": method,
        "related_identifiers": related,
    }
    if meta["zenodo"]["communities"]:
        zenodo["communities"] = [{"identifier": c} for c in meta["zenodo"]["communities"]]
    if meta["funding"]:
        zenodo["grants"] = [
            {"id": f"{g['funder_doi']}::{g['award']}"}
            for g in meta["funding"]
            if g.get("funder_doi") and g.get("award")
        ]
    return zenodo


def dump_cff(cff: dict[str, Any]) -> str:
    """Return ``cff`` as Citation File Format YAML, without YAML anchors or aliases."""
    import yaml

    class _Plain(yaml.SafeDumper):
        def ignore_aliases(self, data: Any) -> bool:
            return True

    return yaml.dump(cff, Dumper=_Plain, sort_keys=False, allow_unicode=True)
