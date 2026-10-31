"""Apply ``mutations:`` and ``modifications:`` blocks, and write them back.

The blocks are plain data, as a config holds them::

    mutations:
      - {site: A:S1500, to: CYS}
      - {chains: [A, C, E], residues: [T301], to: ASN}
    modifications:
      - deprotonate: {site: A:K63, atom_name: NZ}
      - attach:
          type: glycan
          reaction: N-glycosylation
          atom_name: ND2
          fragment_atom_name: C1
          fragment_resnum: 2
          chains: [A, C, E]
          glycan_library: structures/glycans
          glycans: {G55220VL: [N58, N107], G42666HT: [N103]}

A group entry (``chains`` with ``residues`` or ``glycans``, ``sites``, or
``sites_file``) expands to one entry per site, so every site is its own
record. Mutations are applied first, then the modifications in order.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import yaml

from polyzymd.modifying_proteins import api
from polyzymd.modifying_proteins.sites import expand_sites, parse_site

_GROUP_KEYS = ("chains", "residues", "sites", "sites_file", "glycans", "glycan_library")


def expand(spec: dict, base_dir: str | Path = ".") -> tuple[list[dict], list[dict]]:
    """Return the mutations and modifications of ``spec``, one entry per site."""
    base_dir = Path(base_dir)
    mutations = []
    for entry in spec.get("mutations") or []:
        entry = dict(entry)
        target = entry.get("to")
        if isinstance(target, dict) and "ccd_file" in target:
            entry["to"] = {**target, "ccd_file": str(base_dir / target["ccd_file"])}
        for site, _ in _sites(entry, base_dir):
            mutations.append({**_without_groups(entry), "site": site})
    modifications = []
    for entry in spec.get("modifications") or []:
        if not isinstance(entry, dict) or len(entry) != 1:
            raise ValueError(
                f"Each modification is one verb (deprotonate, protonate, attach): {entry!r}"
            )
        ((verb, arguments),) = entry.items()
        if verb == "mutate":
            raise ValueError("Write mutations under `mutations:`; they are applied first.")
        if verb not in ("deprotonate", "protonate", "attach"):
            raise ValueError(
                f"Unknown modification {verb!r}; use deprotonate, protonate or attach."
            )
        arguments = dict(arguments)
        if verb == "attach" and "glycans" in arguments:
            modifications.extend(_glycan_group(arguments, base_dir))
            continue
        sites = _sites(arguments, base_dir)
        if len(sites) > 1 and arguments.get("name"):
            raise ValueError("A group entry cannot set `name`; each site is named in order.")
        for site, fragment in sites:
            single = {**_without_groups(arguments), "site": site}
            if fragment is not None:
                single["fragment"] = fragment
            if verb == "attach":
                single["fragment"] = _resolved_fragment(single, base_dir)
            modifications.append({verb: single})
    return mutations, modifications


def apply_spec(
    protein: api.ModifiedProtein,
    spec: dict,
    base_dir: str | Path = ".",
    relax: str = "joint",
) -> None:
    """Apply the blocks of ``spec`` to ``protein``.

    ``relax`` is ``"joint"`` (place everything, then one relaxation),
    ``"per_modification"`` (mBuild relaxes after each change) or ``"none"``.
    """
    if relax not in ("joint", "per_modification", "none"):
        raise ValueError(f"relax {relax!r}; use joint, per_modification or none.")
    each = relax == "per_modification"
    mutations, modifications = expand(spec, base_dir)
    for entry in mutations:
        api.mutate(
            protein,
            entry["site"],
            entry["to"],
            stereo=entry.get("stereo"),
            relax=each,
            minimize=each,
        )
    for entry in modifications:
        ((verb, arguments),) = entry.items()
        if verb == "attach":
            arguments = dict(arguments)
            fragment = arguments.pop("fragment")
            api.attach(protein, fragment, relax=each, minimize=each, **arguments)
        else:
            getattr(api, verb)(protein, arguments["site"], arguments["atom_name"])
    if relax == "joint":
        api.relax_all(protein)


def to_spec(protein: api.ModifiedProtein, relative_to: str | Path | None = None) -> dict:
    """Return the blocks that rebuild ``protein``'s modifications, one entry per record.

    File paths are written relative to ``relative_to`` (the folder of the
    config that will hold the blocks) when it is given.
    """
    mutations, modifications = [], []
    for record in protein.records:
        site = str(record.site)
        arguments = _paths_relative(_without_hashes(record.arguments), relative_to)
        if record.kind == "mutate":
            mutations.append({"site": site, **arguments})
        elif record.kind == "attach":
            if arguments.get("fragment") is None:
                raise ValueError(
                    f"Attachment {record.name} used a fragment built outside "
                    "polyzymd.modifying_proteins, so it cannot be written as data. Build it "
                    "with fragment_from_smiles, fragment_from_pdb or fragment_from_ccd."
                )
            entry = {"type": record.type, "name": record.name, "site": site, **arguments}
            modifications.append({"attach": entry})
        else:
            modifications.append({record.kind: {"site": site, **arguments}})
    spec = {}
    if mutations:
        spec["mutations"] = mutations
    if modifications:
        spec["modifications"] = modifications
    return spec


def to_yaml(protein: api.ModifiedProtein, relative_to: str | Path | None = None) -> str:
    """Return ``to_spec(protein, relative_to)`` as YAML, ready for a config."""
    return yaml.safe_dump(to_spec(protein, relative_to), sort_keys=False, default_flow_style=None)


def _sites(arguments: dict, base_dir: Path) -> list[tuple[str, dict | None]]:
    """Return each site of an entry with the fragment its sites_file row names."""
    if "site" in arguments:
        return [(str(parse_site(arguments["site"])), None)]
    if "sites" in arguments:
        return [(str(parse_site(site)), None) for site in arguments["sites"]]
    if "chains" in arguments and "residues" in arguments:
        sites = expand_sites(arguments["chains"], arguments["residues"])
        return [(str(site), None) for site in sites]
    if "sites_file" in arguments:
        return _read_sites_file(arguments, base_dir)
    raise ValueError(
        "Give `site`, `sites`, `chains` with `residues`, or `sites_file`: "
        f"{ {k: v for k, v in arguments.items() if k != 'fragment'} }"
    )


def _read_sites_file(arguments: dict, base_dir: Path) -> list[tuple[str, dict | None]]:
    path = base_dir / arguments["sites_file"]
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        header = [name.strip() for name in reader.fieldnames or []]
        rows = [
            {key.strip(): (value or "").strip() for key, value in row.items() if key}
            for row in reader
        ]
    if "site" not in header:
        raise ValueError(f"{path} needs a `site` column (and `glycan` for glycan groups).")
    if not rows:
        raise ValueError(f"{path} has a header but no sites.")
    library = base_dir / arguments.get("glycan_library", ".")
    return [
        (
            str(parse_site(row["site"])),
            {"pdb_path": str(_glycan_file(library, row["glycan"]))} if row.get("glycan") else None,
        )
        for row in rows
    ]


def _glycan_group(arguments: dict, base_dir: Path) -> list[dict]:
    for key in ("name", "site", "sites", "sites_file", "residues"):
        if key in arguments:
            raise ValueError(
                f"A `glycans` group lists its sites under each accession; remove `{key}`."
            )
    library = base_dir / arguments.get("glycan_library", ".")
    chains = arguments.get("chains")
    base = {k: v for k, v in arguments.items() if k not in _GROUP_KEYS and k != "name"}
    entries = []
    for chain in chains or [None]:
        for accession, residues in arguments["glycans"].items():
            residues = [residues] if isinstance(residues, str) else residues
            fragment = {"pdb_path": str(_glycan_file(library, accession))}
            for residue in residues:
                site = parse_site(residue if chain is None else f"{chain}:{residue}")
                entries.append({"attach": {**base, "site": str(site), "fragment": fragment}})
    return entries


def _glycan_file(library: Path, accession: str) -> Path:
    for candidate in (library / f"glycam_{accession}.pdb", library / f"{accession}.pdb"):
        if candidate.is_file():
            return candidate
    raise ValueError(
        f"No file for glycan {accession} in {library}: expected glycam_{accession}.pdb "
        f"or {accession}.pdb with GLYCAM06 residue names and CONECT records."
    )


def _resolved_fragment(entry: dict, base_dir: Path) -> dict:
    fragment = entry.get("fragment")
    if fragment is None:
        raise ValueError(
            f"attach at {entry['site']} has no `fragment` (give smiles + resname, pdb_path, "
            "or ccd + link_atom; a glycan group uses `glycans`)."
        )
    if isinstance(fragment, dict) and "pdb_path" in fragment:
        return {**fragment, "pdb_path": str(base_dir / fragment["pdb_path"])}
    return fragment


def _paths_relative(arguments: dict, relative_to: str | Path | None) -> dict:
    if relative_to is None:
        return arguments
    arguments = dict(arguments)
    for key, field in (("fragment", "pdb_path"), ("fragment", "ccd_file"), ("to", "ccd_file")):
        value = arguments.get(key)
        if isinstance(value, dict) and field in value:
            path = os.path.relpath(Path(value[field]).resolve(), Path(relative_to).resolve())
            arguments[key] = {**value, field: path}
    return arguments


def _without_groups(arguments: dict) -> dict:
    return {k: v for k, v in arguments.items() if k not in _GROUP_KEYS}


def _without_hashes(value):
    if isinstance(value, dict):
        return {k: _without_hashes(v) for k, v in value.items() if k != "sha256"}
    return value
