"""Write src/polyzymd/data/openff_release_dois.json from Zenodo.

"Methods to cite" prints the Zenodo DOI of the OpenFF Toolkit, Interchange and
openff-forcefields release a build used. This script asks Zenodo for every
version of each concept record and writes the version -> DOI table PolyzyMD
ships. Run it before a release (it needs the network):

    python scripts/update_openff_release_dois.py
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

#: Package (as build_provenance.json names it) -> Zenodo concept record.
CONCEPTS = {"toolkit": "597754", "interchange": "8147764", "forcefields": "3473553"}
TABLE = (
    Path(__file__).resolve().parents[1] / "src" / "polyzymd" / "data" / "openff_release_dois.json"
)


def _get(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:  # noqa: S310 - fixed https URL
        return json.load(response)


def versions(concept: str) -> dict[str, str]:
    """Version (without a leading 'v') -> DOI of every release of a Zenodo concept record."""
    latest = _get(f"https://zenodo.org/api/records/{concept}")["id"]
    found: dict[str, str] = {}
    page = 1
    while True:
        hits = _get(f"https://zenodo.org/api/records/{latest}/versions?size=25&page={page}")["hits"]
        for hit in hits["hits"]:
            version = str(hit["metadata"].get("version") or "").removeprefix("v")
            if version:
                found[version] = hit["doi"]
        if not hits["hits"] or page * 25 >= hits["total"]:
            return dict(sorted(found.items()))
        page += 1


def main() -> None:
    table = {
        "source": "Zenodo versions of the concept records "
        + ", ".join(f"10.5281/zenodo.{c}" for c in CONCEPTS.values()),
        **{package: versions(concept) for package, concept in CONCEPTS.items()},
    }
    TABLE.write_text(json.dumps(table, indent=1) + "\n", encoding="utf-8")
    print(f"{TABLE}: " + ", ".join(f"{p} {len(table[p])}" for p in CONCEPTS))


if __name__ == "__main__":
    main()
