"""The ``mutations:``, ``modifications:`` and ``modification_parameters:`` keys of a protein.

A protein solute lists what to change and in which order::

    solutes:
      - type: protein
        pdb_path: structures/trimer.pdb
        force_field: [ff14sb_off_impropers_0.0.4.offxml]
        charges: {method: library}
        mutations:                                # applied first, in order
          - {site: A:S1381, to: 4II}
          - {chains: [A, C, E], residues: [T301], to: ASN}
        modifications:                            # then these, in order
          - deprotonate: {site: A:K63, atom_name: NZ}
          - attach:
              type: glycan
              reaction: N-glycosylation
              chains: [A, C, E]
              glycans: {G55220VL: [N58, N107]}
        modification_parameters:                  # one set per site type
          glycan: {force_field: [GLYCAM_06j-1.xml], charges: {method: library}}

Every group entry (``chains`` with ``residues`` or ``glycans``, ``sites``,
``sites_file``) stands for one modification per site. The models here check
what can be checked without the structure; residue codes are checked against
the PDB when the build plans the modifications.
"""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path
from typing import Any, Literal, Optional, Union

from pydantic import Field, PrivateAttr, field_validator, model_validator

from polyzymd.config.schema import ChargeSpec, _ConfigModel

#: The force field file of GLYCAM06j glycans (OpenMM's amber14 data folder).
GLYCAM_FORCE_FIELD = "GLYCAM_06j-1.xml"

#: Site types that take a parameter set, in the order the docs list them.
SITE_TYPES = ("glycan", "polymer", "generic", "noncanonical")

#: Codes of the 20 canonical residues and their D forms: a mutation to one stays ff14SB.
_CANONICAL_CODES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
    + "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
)
_ONE_LETTER = frozenset("ARNDCQEGHILKMFPSTWYV")


def _check_site(text: str) -> str:
    from polyzymd.modifying_proteins.sites import parse_site

    return str(parse_site(text))


class _Sites(_ConfigModel):
    """Where an entry applies: one ``site``, several ``sites``, ``chains`` × ``residues``,
    or the rows of a ``sites_file`` (CSV with a ``site`` column)."""

    site: Optional[str] = Field(None, description="One residue: A:N58, A:ASN58 or A:N58A")
    sites: Optional[list[str]] = Field(None, min_length=1, description="Several residues")
    chains: Optional[list[str]] = Field(
        None, min_length=1, description="Chains each of `residues` (or `glycans`) applies to"
    )
    residues: Optional[list[str]] = Field(
        None, min_length=1, description="Residues without a chain (N58), with `chains`"
    )
    sites_file: Optional[Path] = Field(
        None, description="CSV with a `site` column (and `glycan` for glycans)"
    )

    @field_validator("site")
    @classmethod
    def _one(cls, value: Optional[str]) -> Optional[str]:
        return None if value is None else _check_site(value)

    @field_validator("sites")
    @classmethod
    def _many(cls, value: Optional[list[str]]) -> Optional[list[str]]:
        return None if value is None else [_check_site(text) for text in value]

    def _site_forms(self) -> list[str]:
        return [
            form
            for form, given in (
                ("site", self.site is not None),
                ("sites", self.sites is not None),
                ("chains + residues", self.residues is not None),
                ("sites_file", self.sites_file is not None),
            )
            if given
        ]

    def _require_one_site_form(self) -> None:
        forms = self._site_forms()
        if len(forms) != 1:
            raise ValueError(
                "give exactly one of `site`, `sites`, `chains` with `residues`, or "
                f"`sites_file`; got {', '.join(forms) or 'none'}"
            )
        if self.residues is not None and self.chains is None:
            raise ValueError("`residues` needs `chains` (or write full sites: A:N58)")
        if self.chains is not None and self.residues is None:
            raise ValueError("`chains` needs `residues` (the residues to change on each chain)")


class SmilesSideChain(_ConfigModel):
    """A side chain from SMILES; ``*`` marks the atom that bonds to CA."""

    smiles: str = Field(..., description="Side chain SMILES with one * (the bond to CA)")
    resname: str = Field(
        ..., min_length=1, max_length=3, description="Residue name, 1-3 characters"
    )


class CcdFile(_ConfigModel):
    """A residue from a local CCD definition (.cif), for codes mBuild does not bundle."""

    ccd_file: Path = Field(..., description="Path to the component's CCD .cif file")


class Mutation(_Sites):
    """One entry of ``mutations:``: replace the side chain at each site."""

    to: Union[str, SmilesSideChain, CcdFile] = Field(
        ..., description="A CCD code (ALA, 4II, DLE), {smiles, resname}, or {ccd_file}"
    )
    stereo: Optional[Literal["L", "D"]] = Field(
        None, description="Handedness after the mutation (default: the component's own)"
    )

    @model_validator(mode="after")
    def _where(self) -> "Mutation":
        self._require_one_site_form()
        return self

    @property
    def is_canonical(self) -> bool:
        """Whether the mutation keeps a canonical residue (ff14SB, no site)."""
        return isinstance(self.to, str) and self.to.upper() in _CANONICAL_CODES


class Protonation(_Sites):
    """Arguments of ``protonate`` or ``deprotonate``: the atom that gains or loses a proton."""

    atom_name: str = Field(..., description="The atom that gains or loses H (NZ, OD2, N)")

    @model_validator(mode="after")
    def _where(self) -> "Protonation":
        self._require_one_site_form()
        return self


class GraftPolymer(_ConfigModel):
    """A polymer chain built for an attachment of ``type: polymer``.

    The monomers and caps come from a ``type: polymer`` solute
    (``monomers_from: <its name>``) or are written here as for a solute.
    ``length`` counts the two caps, and positions (in ``sequence`` and the
    attachment's ``attach_at``) count the whole chain, caps included. The
    chain is drawn anew for each replicate; ``random_seed`` fixes one chain
    for every replicate.
    """

    monomers_from: Optional[str] = Field(None, description="Name of a type: polymer solute")
    monomers: Optional[list[dict[str, Any]]] = Field(None, description="Monomers, as for a solute")
    end_groups: Optional[dict[str, Any]] = Field(None, description="Caps, as for a solute")
    length: int = Field(..., ge=3, description="Units of the chain, the two caps included")
    sequence: Optional[dict[str, Any]] = Field(None, description="As for a solute")
    random_seed: Optional[int] = Field(
        None, description="One chain for every replicate (default: drawn per replicate)"
    )
    _resolved: Optional[dict[str, Any]] = PrivateAttr(None)

    @model_validator(mode="after")
    def _monomers_given_once(self) -> "GraftPolymer":
        if (self.monomers_from is None) == (self.monomers is None):
            raise ValueError(
                "a polymer fragment takes its monomers from a polymer solute "
                "(monomers_from: <name>) or lists them (monomers + end_groups), not both"
            )
        if self.monomers is not None:
            if self.end_groups is None:
                raise ValueError("a polymer fragment that lists monomers needs end_groups")
            self.resolve(self.monomers, self.end_groups)
        elif self.end_groups is not None:
            raise ValueError("with monomers_from the caps come from that solute too")
        return self

    def resolve(self, monomers: list[dict[str, Any]], end_groups: dict[str, Any]) -> None:
        """Check and keep the chain's monomers and caps (from this entry or a solute)."""
        from polyzymd.config.schema import MonomerSpec, PolymerEndGroups, PolymerSequence

        specs = [MonomerSpec.model_validate(m) for m in monomers]
        missing = [m.label for m in specs if not m.smiles or not m.name]
        if missing:
            raise ValueError(
                f"grafted chain monomers {missing} need `name` and star-sited `smiles`"
            )
        caps = PolymerEndGroups.model_validate(end_groups)
        label_of = {m.label: m.label for m in specs}
        label_of.update({m.name: m.label for m in specs})
        cap_labels = {
            end: {
                **{o.label: o.label for o in getattr(caps, end)},
                **{o.name: o.label for o in getattr(caps, end) if o.name},
            }
            for end in ("head", "tail")
        }
        PolymerSequence.model_validate(self.sequence or {}).pattern(
            self.length, label_of, cap_labels
        )
        self._resolved = {
            "monomers": [m.model_dump(exclude_none=True) for m in specs],
            "end_groups": caps.model_dump(exclude_none=True),
        }

    def as_source(self, force_field: list[str]) -> dict[str, Any]:
        """Return the chain as ``polyzymd.modifying_proteins`` takes it."""
        if self._resolved is None:
            raise ValueError(
                f"monomers_from {self.monomers_from!r} was not resolved; it names a polymer "
                "solute of the same config"
            )
        source: dict[str, Any] = {**self._resolved, "length": self.length}
        if self.sequence:
            source["sequence"] = self.sequence
        if self.random_seed is not None:
            source["random_seed"] = self.random_seed
        source["force_field"] = list(force_field)
        return source


class Fragment(_ConfigModel):
    """The group an attachment bonds: SMILES, a PDB file, a CCD component or a polymer."""

    smiles: Optional[str] = Field(None, description="SMILES with * at the bond site")
    resname: Optional[str] = Field(None, max_length=3, description="Residue name with smiles")
    pdb_path: Optional[Path] = Field(None, description="PDB file with CONECT records")
    charge: int = Field(0, description="Net formal charge of a pdb_path fragment")
    ccd: Optional[str] = Field(None, description="CCD code of the fragment")
    ccd_file: Optional[Path] = Field(
        None, description="Local CCD definition (.cif) of the fragment, named <CODE>.cif"
    )
    link_atom: Optional[str] = Field(None, description="Atom of the CCD fragment that bonds")
    polymer: Optional[GraftPolymer] = Field(None, description="A chain built for a polymer graft")

    @model_validator(mode="after")
    def _one_source(self) -> "Fragment":
        given = [
            name
            for name, present in (
                ("smiles", self.smiles is not None),
                ("pdb_path", self.pdb_path is not None),
                ("ccd", self.ccd is not None),
                ("ccd_file", self.ccd_file is not None),
                ("polymer", self.polymer is not None),
            )
            if present
        ]
        if len(given) != 1:
            raise ValueError(
                "a fragment is one of {smiles, resname}, {pdb_path}, {ccd, link_atom}, "
                f"{{ccd_file, link_atom}} or {{polymer}}; got {', '.join(given) or 'none'}"
            )
        if self.smiles is not None and not self.resname:
            raise ValueError("a SMILES fragment needs `resname` (1-3 characters)")
        if (self.ccd is not None or self.ccd_file is not None) and not self.link_atom:
            raise ValueError("a CCD fragment needs `link_atom`, the atom that bonds")
        return self

    def as_source(self, polymer_force_field: Optional[list[str]] = None) -> dict[str, Any]:
        """Return the fragment as polyzymd.modifying_proteins takes it."""
        if self.polymer is not None:
            return {"polymer": self.polymer.as_source(list(polymer_force_field or []))}
        if self.smiles is not None:
            return {"smiles": self.smiles, "resname": self.resname}
        if self.pdb_path is not None:
            source: dict[str, Any] = {"pdb_path": str(self.pdb_path)}
            if self.charge:
                source["charge"] = self.charge
            return source
        if self.ccd_file is not None:
            return {"ccd_file": str(self.ccd_file), "link_atom": self.link_atom}
        return {"ccd": self.ccd, "link_atom": self.link_atom}


class Attachment(_Sites):
    """Arguments of ``attach``: bond a fragment to an atom of each site's residue.

    ``type`` (glycan, polymer or generic) is never inferred: it picks the
    segment the new atoms get, the parameter set that covers them and the
    checks that apply. ``reaction`` (an mBuild named reaction or a reaction
    SMARTS) decides which atoms leave and which bonds form; without it, the
    atom-level fields do.
    """

    type: Literal["glycan", "polymer", "generic"] = Field(
        ..., description="glycan, polymer or generic; picks segment, parameters and checks"
    )
    name: Optional[str] = Field(
        None, max_length=4, description="Segment name (default G001, P001, M001 in order)"
    )
    atom_name: Optional[str] = Field(
        None,
        description="Atom of the protein residue that bonds (ND2, NZ); a glycan defaults to "
        "Asn ND2, Ser OG or Thr OG1",
    )
    fragment: Optional[Fragment] = Field(None, description="The group to bond")
    glycans: Optional[dict[str, Union[list[str], str]]] = Field(
        None, description="Glycan accession: residues it goes on ({G42666HT: [N103]})"
    )
    glycan_library: Optional[Path] = Field(
        None, description="Folder of glycam_<accession>.pdb files (default structures/glycans)"
    )
    reaction: Optional[str] = Field(None, description="mBuild named reaction or reaction SMARTS")
    fragment_atom_name: Optional[str] = Field(
        None, description="Atom of the fragment that bonds (a glycan defaults to C1)"
    )
    fragment_resnum: Optional[int] = Field(
        None,
        description="Residue of the fragment that bonds (a glycan defaults to the first "
        "residue of its file with a C1)",
    )
    fragment_resname: Optional[str] = Field(None, description="Residue name of the fragment")
    leaving_atom_names: Optional[Union[str, list[str]]] = Field(
        None, description="Protein atoms that leave"
    )
    fragment_leaving_atom_names: Optional[Union[str, list[str]]] = Field(
        None, description="Fragment atoms that leave"
    )
    merge: bool = Field(False, description="Merge the fragment into the residue (ring closures)")
    bond_order: int = Field(1, ge=1, le=3, description="Order of the new bond")
    attach_at: Optional[int] = Field(
        None,
        ge=1,
        description="Position of a polymer fragment that bonds, 1..length (caps 1, length)",
    )

    @model_validator(mode="before")
    @classmethod
    def _no_per_entry_parameters(cls, data: Any) -> Any:
        if isinstance(data, dict) and "parameters" in data:
            raise ValueError(
                "one parameter set per attachment type in 1.4: set it under "
                "modification_parameters (glycan, polymer, generic, noncanonical) and "
                f"remove `parameters` from this {data.get('type', '')} entry"
            )
        return data

    @field_validator("glycans")
    @classmethod
    def _glycan_residues(cls, value):
        if value is None:
            return value
        return {
            accession: [residues] if isinstance(residues, str) else list(residues)
            for accession, residues in value.items()
        }

    @model_validator(mode="after")
    def _where_and_what(self) -> "Attachment":
        if self.atom_name is None and self.type != "glycan":
            raise ValueError(
                f"a {self.type} attachment needs atom_name, the protein atom that bonds; only "
                "a glycan's defaults (Asn ND2, Ser OG, Thr OG1)"
            )
        if self.glycans is not None:
            if self.type != "glycan":
                raise ValueError(f"`glycans` is for type: glycan, not type: {self.type}")
            if self.site is not None or self.sites is not None or self.sites_file is not None:
                raise ValueError(
                    "a `glycans` group names its residues itself; remove site/sites/sites_file"
                )
            if self.residues is not None:
                raise ValueError("a `glycans` group lists its residues under each accession")
            if self.fragment is not None:
                raise ValueError("a `glycans` group takes its fragments from glycan_library")
        else:
            self._require_one_site_form()
            if self.fragment is None and not (self.sites_file and self.type == "glycan"):
                raise ValueError(
                    "an attachment needs `fragment` ({smiles, resname}, {pdb_path} or "
                    "{ccd, link_atom}); a glycan group may use `glycans:` or a sites_file "
                    "with a `glycan` column instead"
                )
        polymer = self.fragment is not None and self.fragment.polymer is not None
        if polymer and self.type != "polymer":
            raise ValueError(f"a polymer fragment is attached with type: polymer, not {self.type}")
        if polymer and self.attach_at is None:
            raise ValueError(
                "a grafted chain needs `attach_at`: the position that bonds, 1..length with the "
                "caps at 1 and length"
            )
        if self.attach_at is not None and not polymer:
            raise ValueError("`attach_at` is for polymer fragments; use fragment_resnum otherwise")
        if polymer and self.fragment_resnum is not None:
            raise ValueError("give attach_at for a grafted chain, not fragment_resnum")
        if polymer and self.attach_at is not None:
            length = self.fragment.polymer.length
            if self.attach_at > length:
                raise ValueError(
                    f"attach_at {self.attach_at} is outside the chain's positions 1..{length} "
                    f"(the caps are 1 and {length})"
                )
        if self.name is not None and self._is_group:
            raise ValueError("a group entry cannot set `name`; its sites are named in order")
        return self

    @property
    def _is_group(self) -> bool:
        return self.site is None


class Modification(_ConfigModel):
    """One entry of ``modifications:``: a single verb and its arguments."""

    deprotonate: Optional[Protonation] = None
    protonate: Optional[Protonation] = None
    attach: Optional[Attachment] = None

    @model_validator(mode="before")
    @classmethod
    def _one_verb(cls, data: Any) -> Any:
        if isinstance(data, dict):
            if "mutate" in data:
                raise ValueError(
                    "mutations go under `mutations:` (they are applied before every "
                    "modification): move this entry there as {site: ..., to: ...}"
                )
            unknown = sorted(set(data) - {"deprotonate", "protonate", "attach"})
            if unknown:
                raise ValueError(
                    f"unknown modification {', '.join(unknown)}; use deprotonate, protonate "
                    "or attach"
                )
            if len(data) != 1:
                raise ValueError(
                    "each modification is one verb (deprotonate, protonate or attach); "
                    "write one per entry"
                )
        return data

    @property
    def verb(self) -> str:
        return next(name for name in ("deprotonate", "protonate", "attach") if getattr(self, name))

    @property
    def arguments(self) -> Union[Protonation, Attachment]:
        return getattr(self, self.verb)


class ModificationParameterSet(_ConfigModel):
    """The force field and charges of one site type.

    An OpenFF force field (e.g. openff-2.3.0.offxml) takes NAGL charges on a
    capped model of each site. GLYCAM06j takes its own library charges.
    """

    force_field: list[str] = Field(..., min_length=1, description="Force field files")
    charges: ChargeSpec = Field(..., description="nagl (OpenFF force fields) or library (GLYCAM)")

    @model_validator(mode="after")
    def _consistent(self) -> "ModificationParameterSet":
        if self.is_glycam:
            if self.force_field != [GLYCAM_FORCE_FIELD] or self.charges.method != "library":
                raise ValueError(
                    f"GLYCAM06j is force_field: [{GLYCAM_FORCE_FIELD}] with charges: "
                    "{method: library}; its residue templates carry both"
                )
        elif self.charges.method != "nagl":
            raise ValueError(
                "site charges come from NAGL on a capped model of each site: use "
                "charges: {method: nagl, model: <model .pt file>}"
            )
        return self

    @property
    def is_glycam(self) -> bool:
        return any(Path(name).name == GLYCAM_FORCE_FIELD for name in self.force_field)


class ModificationParameters(_ConfigModel):
    """The parameter set of each site type; a type with a site must have one."""

    glycan: Optional[ModificationParameterSet] = Field(None, description="Glycan sites")
    polymer: Optional[ModificationParameterSet] = Field(None, description="Grafted polymers")
    generic: Optional[ModificationParameterSet] = Field(None, description="Dyes, labels, clicks")
    noncanonical: Optional[ModificationParameterSet] = Field(
        None, description="Non-canonical residues (mutations, protonation states)"
    )

    @model_validator(mode="after")
    def _glycam_for_glycans_only(self) -> "ModificationParameters":
        misused = [
            name
            for name in ("polymer", "generic", "noncanonical")
            if getattr(self, name) is not None and getattr(self, name).is_glycam
        ]
        if misused:
            raise ValueError(
                f"GLYCAM06j parameterizes glycans only; {', '.join(misused)} cannot use it"
            )
        return self

    def as_dict(self) -> dict[str, ModificationParameterSet]:
        return {name: getattr(self, name) for name in SITE_TYPES if getattr(self, name)}


def required_site_types(
    mutations: list[Mutation], modifications: list[Modification]
) -> list[tuple[str, str]]:
    """Return (site type, why) for each type the entries are known to create sites of.

    A protonation change outside the ff14SB variants also makes a
    non-canonical site; that is found when the build plans the change.
    """
    needed = []
    for mutation in mutations:
        if not mutation.is_canonical:
            needed.append(("noncanonical", f"the mutation to {_describe_target(mutation.to)}"))
    for modification in modifications:
        if modification.verb == "attach":
            needed.append((modification.attach.type, f"a {modification.attach.type} attachment"))
    return needed


def _describe_target(target) -> str:
    if isinstance(target, str):
        return target
    if isinstance(target, SmilesSideChain):
        return f"the SMILES side chain {target.resname}"
    return str(target.ccd_file)


def declared_modifications(
    mutations: list[Mutation],
    modifications: list[Modification],
    parameters: Optional[ModificationParameters],
) -> list[dict[str, Any]]:
    """Return what the config asks for, one entry per site, read from the config alone.

    This is the declared side the audit compares the built protein with, so
    it expands the group entries itself instead of reusing the builder's
    expansion: a wrong expansion in the build then shows as a mismatch.
    """
    sets = parameters.as_dict() if parameters else {}
    declared = []
    for n, mutation in enumerate(mutations):
        for site in _expand(mutation):
            entry = {
                "entry": f"mutations.{n}",
                "kind": "mutate",
                "site": site,
                "to": _target_data(mutation.to),
            }
            if mutation.stereo:
                entry["stereo"] = mutation.stereo
            if not mutation.is_canonical:
                entry["parameters"] = _set_data(sets.get("noncanonical"))
            declared.append(entry)
    for n, modification in enumerate(modifications):
        verb, arguments = modification.verb, modification.arguments
        where = f"modifications.{n}.{verb}"
        if verb != "attach":
            for site in _expand(arguments):
                declared.append(
                    {"entry": where, "kind": verb, "site": site, "atom_name": arguments.atom_name}
                )
            continue
        # Every field that changes the chemistry: all of the entry but how it names
        # its sites (and the fragment, replaced by its source data and hash).
        fields = {
            key: value
            for key, value in arguments.model_dump(mode="json", exclude=_SITE_FORM_KEYS).items()
            if value is not None
        }
        for key in ("leaving_atom_names", "fragment_leaving_atom_names"):
            if isinstance(fields.get(key), str):
                fields[key] = [fields[key]]
        polymer_set = sets.get("polymer")
        polymer_force_field = list(polymer_set.force_field) if polymer_set else []
        for site, fragment in _attachment_sites(arguments, polymer_force_field):
            declared.append(
                {
                    "entry": where,
                    "kind": "attach",
                    "site": site,
                    **fields,
                    "fragment": fragment,
                    "parameters": _set_data(sets.get(arguments.type)),
                }
            )
    return declared


#: Keys that say where an entry applies, not what it does.
_SITE_FORM_KEYS = {
    "site",
    "sites",
    "chains",
    "residues",
    "sites_file",
    "glycans",
    "glycan_library",
    "name",
    "fragment",
}


def modification_fingerprint(solute: Any) -> dict[str, Any]:
    """Everything about a protein's modifications that changes the system, as plain data.

    One entry per site with every chemistry field (files by SHA-256, not path),
    the parameter sets, ``max_charge_shift`` and ``relax``; empty for an
    unmodified solute. The analysis identity and the run-folder token hash it.
    """
    if not getattr(solute, "is_modified", False):
        return {}

    def portable(value: Any) -> Any:
        if isinstance(value, dict):
            return {
                k: portable(v)
                for k, v in value.items()
                if k not in ("pdb_path", "ccd_file", "entry")
            }
        if isinstance(value, list):
            return [portable(v) for v in value]
        return value

    declared = declared_modifications(
        solute.mutations, solute.modifications, solute.modification_parameters
    )
    return {
        "modifications": [portable(item) for item in declared],
        "max_charge_shift": solute.max_charge_shift,
        "relax": solute.relax,
    }


def modification_identity(solute: Any) -> str:
    """SHA-256 (hex) of :func:`modification_fingerprint`; "" for an unmodified solute."""
    import json

    fingerprint = modification_fingerprint(solute)
    if not fingerprint:
        return ""
    encoded = json.dumps(fingerprint, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _expand(entry: _Sites) -> list[str]:
    from polyzymd.modifying_proteins.sites import parse_site

    if entry.site is not None:
        return [entry.site]
    if entry.sites is not None:
        return list(entry.sites)
    if entry.residues is not None:
        return [str(parse_site(f"{c}:{r}")) for c in entry.chains for r in entry.residues]
    return [str(parse_site(row["site"])) for row in _csv_rows(entry.sites_file)]


def _attachment_sites(
    attachment: Attachment, polymer_force_field: Optional[list[str]] = None
) -> list[tuple[str, dict[str, Any]]]:
    from polyzymd.modifying_proteins.sites import parse_site

    library = attachment.glycan_library
    if attachment.glycans is not None:
        return [
            (
                str(parse_site(residue if chain is None else f"{chain}:{residue}")),
                _file_data(_glycan_path(library, accession)),
            )
            for chain in attachment.chains or [None]
            for accession, residues in attachment.glycans.items()
            for residue in residues
        ]
    if attachment.sites_file is not None and attachment.fragment is None:
        return [
            (str(parse_site(row["site"])), _file_data(_glycan_path(library, _glycan_of(row))))
            for row in _csv_rows(attachment.sites_file)
        ]
    fragment = _fragment_data(attachment.fragment, polymer_force_field)
    return [(site, fragment) for site in _expand(attachment)]


def _glycan_of(row: dict[str, str]) -> str:
    if not row.get("glycan"):
        raise ValueError(
            f"sites file row {row} names no glycan; a glycan group read from a sites file "
            "needs a `glycan` column (an accession)"
        )
    return row["glycan"]


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return [
            {key.strip(): (value or "").strip() for key, value in row.items() if key}
            for row in csv.DictReader(handle)
        ]


def glycan_file(library: Optional[Path], accession: str) -> Path:
    """The file of glycan ``accession`` in ``library``: ``glycam_<ID>.pdb`` or ``<ID>.pdb``
    (the first when neither exists)."""
    return _glycan_path(library, accession)


def _glycan_path(library: Optional[Path], accession: str) -> Path:
    folder = Path(library) if library is not None else Path("structures/glycans")
    for candidate in (folder / f"glycam_{accession}.pdb", folder / f"{accession}.pdb"):
        if candidate.is_file():
            return candidate
    return folder / f"glycam_{accession}.pdb"


def _file_data(path: Path) -> dict[str, Any]:
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"pdb_path": str(path), "sha256": digest}


def _fragment_data(
    fragment: Fragment, polymer_force_field: Optional[list[str]] = None
) -> dict[str, Any]:
    if fragment.pdb_path is not None:
        return _file_data(fragment.pdb_path)
    if fragment.ccd_file is not None:
        digest = _file_data(fragment.ccd_file)["sha256"]
        return {
            "ccd_file": str(fragment.ccd_file),
            "sha256": digest,
            "link_atom": fragment.link_atom,
        }
    return fragment.as_source(polymer_force_field)


def resolve_grafts(solutes: list[Any]) -> None:
    """Give each grafted chain with ``monomers_from`` the monomers and caps of that solute."""
    polymers = {solute.name: solute for solute in solutes if solute.type == "polymer"}
    for solute in solutes:
        for modification in getattr(solute, "modifications", None) or []:
            attachment = modification.attach
            if attachment is None or attachment.fragment is None:
                continue
            chain = attachment.fragment.polymer
            if chain is None or chain.monomers_from is None:
                continue
            source = polymers.get(chain.monomers_from)
            if source is None:
                known = ", ".join(sorted(polymers)) or "none"
                raise ValueError(
                    f"monomers_from {chain.monomers_from!r} names no polymer solute "
                    f"(polymer solutes: {known})"
                )
            if source.end_groups is None or any(m.smiles is None for m in source.monomers):
                raise ValueError(
                    f"monomers_from {chain.monomers_from!r}: a grafted chain is built from "
                    "star-sited SMILES and end_groups, which that solute does not give "
                    "(use generation_mode: dynamic there, or list the monomers here)"
                )
            chain.resolve(
                [monomer.model_dump(exclude_none=True) for monomer in source.monomers],
                source.end_groups.model_dump(exclude_none=True),
            )


def _target_data(target) -> Any:
    if isinstance(target, str):
        return target.upper()
    if isinstance(target, SmilesSideChain):
        return {"smiles": target.smiles, "resname": target.resname}
    return _file_data(target.ccd_file) | {"ccd_file": str(target.ccd_file)}


def _set_data(parameter_set: Optional[ModificationParameterSet]) -> Optional[dict[str, Any]]:
    if parameter_set is None:
        return None
    return {
        "force_field": list(parameter_set.force_field),
        "charges": parameter_set.charges.model_dump(exclude_none=True),
    }


def declared_for_audit(config: Any) -> Optional[dict[str, dict[str, Any]]]:
    """What each modified protein of ``config`` declares, for ``modifications.declared``.

    ``{solute name: {"spec", "parameter_sets", "base_dir"}}``, or None when no
    protein has modifications.
    """
    if config is None:
        return None
    declared = {}
    for solute in getattr(config, "solutes", None) or []:
        if not getattr(solute, "is_modified", False):
            continue
        sets = solute.modification_parameters.as_dict() if solute.modification_parameters else {}
        declared[solute.name] = {
            # The config's own expansion, one entry per site; not the builder's.
            "entries": declared_modifications(
                solute.mutations, solute.modifications, solute.modification_parameters
            ),
            "parameter_sets": {name: _set_data(chosen) for name, chosen in sets.items()},
            "base_dir": str(getattr(config, "_config_dir", None) or "."),
        }
    return declared or None


_RESIDUE_NUMBER = re.compile(r"\b(resid|resnum|resSeq|residue)\b")
_SEGMENT_SCOPE = re.compile(r"\b(segid|protein)\b")


def has_attachments(solute: Any) -> bool:
    """Whether a protein solute attaches anything (so it has segments beside PROT)."""
    return any(entry.verb == "attach" for entry in getattr(solute, "modifications", None) or [])


def unscoped_residue_selection(selection: str) -> bool:
    """Whether ``selection`` picks residues by number without naming a segment.

    A modified protein's attachments are segments of their own whose residues
    are numbered from 1, so ``resid 1-10`` also matches glycan, graft and
    label residues; ``segid PROT and resid 1-10`` (or ``protein and ...``) does not.
    """
    return bool(_RESIDUE_NUMBER.search(selection)) and not _SEGMENT_SCOPE.search(selection)
