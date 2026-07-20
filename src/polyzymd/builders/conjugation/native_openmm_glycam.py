"""Native OpenMM GLYCAM handoff for strict N-linked glycan builds.

The exact bundle exists because OpenFF Interchange cannot yet preserve the
complete OpenMM ``NonbondedForce`` exception table needed by GLYCAM. GLYCAM uses
explicit exception parameters whose global 1-4 semantics are lossy when coerced
through ordinary Interchange scaling. PolyzyMD therefore keeps native OpenMM as
the authoritative system and carries an exact sidecar for export until upstream
Interchange can import and export explicit OpenMM exceptions directly.

Covalently attached Sage polymers across Amber/GLYCAM boundaries remain
unsupported. Such systems lack a single force-field owner for cross-boundary
bonded terms and exceptions, so strict routing fails closed instead of inventing
provenance.
"""

from __future__ import annotations

import json
import logging
import re
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.exporters.exact_openmm import create_exact_export_bundle

LOGGER = logging.getLogger(__name__)

GLYCAM_XML_STACK: tuple[str, ...] = (
    "amber14/protein.ff14SB.xml",
    "amber14/GLYCAM_06j-1.xml",
    "amber14/tip3p.xml",
)
NATIVE_GLYCAM_AUDIT_NAME = "native_openmm_glycam_audit.json"
NATIVE_GLYCAM_NONBONDED_CUTOFF_NM = 1.0
NATIVE_GLYCAM_CONSTRAINTS = "HBonds"
NATIVE_GLYCAM_RIGID_WATER = True
NATIVE_GLYCAM_ROUTE_INVARIANTS = {
    "nonbonded_cutoff_nm": NATIVE_GLYCAM_NONBONDED_CUTOFF_NM,
    "constraints": NATIVE_GLYCAM_CONSTRAINTS,
    "rigid_water": NATIVE_GLYCAM_RIGID_WATER,
}
SUPPORTED_WATER_RESIDUES = {"HOH", "WAT", "H2O"}
SUPPORTED_ION_RESIDUES = {"NA", "CL", "K", "MG", "CA"}


def native_glycam_enabled(config: Any) -> bool:
    """Return whether a config explicitly requests native OpenMM GLYCAM."""
    force_field = getattr(config, "force_field", None)
    return (
        getattr(force_field, "conjugate_parameterization", None) == "native_openmm_glycam"
        and getattr(force_field, "glycan_policy", "sage_fallback") == "strict_glycam"
    )


def create_native_openmm_glycam_handoff(
    builder: Any,
    *,
    config: Any,
    construction: Any,
    output_dir: Path | str,
) -> Any:
    """Create a direct OpenMM GLYCAM handoff for a solvated N-linked glycan.

    Parameters
    ----------
    builder : Any
        Prepared system builder containing the solvated authoritative OpenFF
        topology from Pablo and PolyzyMD solvation.
    config : Any
        Simulation configuration that explicitly opted into native GLYCAM.
    construction : Any
        Conjugate construction result with attachment/crosslink provenance.
    output_dir : pathlib.Path or str
        Directory where the audit JSON will be written.

    Returns
    -------
    Any
        Exact OpenMM export bundle with topology, system, positions, sidecar,
        private baseline Interchange, and audit path.
    """
    from openmm import unit
    from openmm.app import PME, HBonds

    _validate_native_glycam_mvp_config(config)
    solvated_topology = getattr(builder, "_solvated_topology", None)
    if solvated_topology is None:
        raise RuntimeError("Native GLYCAM handoff requires a solvated topology")
    _annotate_force_field_domains_from_config(solvated_topology, construction, config)

    sage_residue_names = _sage_residue_names_from_config(config)
    converted = _openff_topology_to_openmm_for_glycam(
        solvated_topology,
        construction,
        sage_residue_names=sage_residue_names,
    )
    force_field = _load_native_glycam_force_field()
    sage_components = _register_disconnected_sage_template_generator(
        force_field, solvated_topology, config, construction
    )
    template_matches = _preflight_glycam_template_matches(force_field, converted.topology)

    asx_residues = _converted_asx_residues(converted)
    crosslink_pairs = _converted_crosslink_pairs(converted)
    system = force_field.createSystem(
        converted.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=NATIVE_GLYCAM_NONBONDED_CUTOFF_NM * unit.nanometer,
        constraints=HBonds,
        rigidWater=NATIVE_GLYCAM_RIGID_WATER,
        residueTemplates=dict.fromkeys(asx_residues, "NLN"),
    )
    audit = _build_native_glycam_audit(
        converted.topology,
        system,
        asx_residues,
        crosslink_pairs,
        renamed_atoms=converted.renamed_atoms,
        template_matches=template_matches,
        sage_components=sage_components,
        construction=construction,
    )
    _require_essential_linkage_terms(audit)
    audit_path = Path(output_dir) / NATIVE_GLYCAM_AUDIT_NAME
    audit_path.write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    exact_export_bundle = create_exact_export_bundle(
        topology=converted.topology,
        system=system,
        positions=converted.positions,
        output_dir=output_dir,
        audit_path=audit_path,
        audit=audit,
    )
    builder._exact_export_bundle = exact_export_bundle
    return exact_export_bundle


@dataclass(frozen=True)
class _ConvertedTopology:
    """OpenMM topology conversion result for native GLYCAM."""

    topology: Any
    positions: Any
    asx_residues: tuple[Any, ...]
    crosslink_pairs: tuple[tuple[Any, Any], ...]
    renamed_atoms: tuple[dict[str, Any], ...]
    sage_template_units: tuple[dict[str, Any], ...] = ()

    @property
    def asx_residue(self) -> Any:
        """Return the single ASX residue for legacy tests and callers."""
        if len(self.asx_residues) != 1:
            raise ValueError(f"Expected one ASX/NLN residue; found {len(self.asx_residues)}")
        return self.asx_residues[0]

    @property
    def crosslink_atoms(self) -> tuple[Any, Any]:
        """Return the single crosslink for legacy tests and callers."""
        if len(self.crosslink_pairs) != 1:
            raise ValueError(f"Expected one GLYCAM crosslink; found {len(self.crosslink_pairs)}")
        return self.crosslink_pairs[0]


def _converted_asx_residues(converted: Any) -> tuple[Any, ...]:
    """Return ASX/NLN residues from current or legacy converted topology objects."""
    if hasattr(converted, "asx_residues"):
        return tuple(converted.asx_residues)
    return (converted.asx_residue,)


def _converted_crosslink_pairs(converted: Any) -> tuple[tuple[Any, Any], ...]:
    """Return crosslink pairs from current or legacy converted topology objects."""
    if hasattr(converted, "crosslink_pairs"):
        return tuple(converted.crosslink_pairs)
    return (converted.crosslink_atoms,)


def _validate_native_glycam_mvp_config(config: Any) -> None:
    """Reject unsupported components before native GLYCAM parameterization."""
    if getattr(config, "engine", None) not in {"openmm", "gromacs"}:
        raise ValueError("Native GLYCAM mode supports only engine='openmm' or engine='gromacs'")
    polymers = getattr(config, "polymers", None)
    if polymers is not None and getattr(polymers, "enabled", False):
        LOGGER.info("Native GLYCAM mode will admit only disconnected precharged Sage polymers")
    solvent = getattr(config, "solvent", None)
    if solvent is None:
        raise ValueError("Native GLYCAM mode requires explicit TIP3P water solvent")
    primary = getattr(solvent, "primary", None)
    if getattr(primary, "type", "water") != "water" or str(
        getattr(primary, "model", "tip3p")
    ) not in {
        "tip3p",
        "WaterModel.TIP3P",
    }:
        raise ValueError("Native GLYCAM mode supports only TIP3P water")
    ions = getattr(solvent, "ions", None)
    if ions is not None and (
        getattr(ions, "neutralize", True)
        or float(getattr(ions, "nacl_concentration", 0.0)) != 0.0
        or float(getattr(ions, "kcl_concentration", 0.0)) != 0.0
        or float(getattr(ions, "mgcl2_concentration", 0.0)) != 0.0
    ):
        raise ValueError(
            "Native GLYCAM mode supports only ions that match Amber/TIP3P templates; "
            "automatic ion placement and neutralization are not yet audited for this route"
        )


def _load_native_glycam_force_field() -> Any:
    """Load the OpenMM-vendored Amber14 ff14SB, GLYCAM, and TIP3P XML stack."""
    from openmm.app import ForceField

    return ForceField(*GLYCAM_XML_STACK)


def _annotate_force_field_domains_from_config(
    topology: Any, construction: Any, config: Any
) -> None:
    """Copy config-level attachment domains onto OpenFF atom metadata when possible."""
    _ = construction
    conjugation = getattr(config, "conjugation", None)
    sage_residue_names = _sage_residue_names_from_config(config)
    residue_domains: dict[str, str] = {}
    if conjugation is None:
        pass
    else:
        for attachment in getattr(conjugation, "attachments", ()) or ():
            if not getattr(attachment, "enabled", True):
                continue
            domain = getattr(getattr(attachment, "moiety", None), "force_field_domain", None)
            if domain is None:
                continue
            product_residues = getattr(
                getattr(attachment, "mechanism", None), "product_residues", None
            )
            for residue_name in (
                getattr(product_residues, "site", None),
                getattr(product_residues, "moiety", None),
                getattr(getattr(attachment, "moiety", None), "residue_name", None),
            ):
                if residue_name:
                    residue_domains[str(residue_name).strip().upper()] = str(domain).strip().lower()
    residue_domains.update(dict.fromkeys(sage_residue_names, "sage"))
    if not residue_domains:
        return
    for molecule in getattr(topology, "molecules", ()):
        molecule_domains = set()
        molecule_residue_names = set()
        for atom in getattr(molecule, "atoms", ()):
            metadata = getattr(atom, "metadata", None)
            if not isinstance(metadata, dict):
                continue
            residue_name = str(metadata.get("residue_name", "")).strip().upper()
            if residue_name:
                molecule_residue_names.add(residue_name[:3])
            domain = residue_domains.get(residue_name) or residue_domains.get(residue_name[:3])
            if domain is None:
                continue
            metadata.setdefault("force_field_domain", domain)
            molecule_domains.add(domain)
        if molecule_residue_names & sage_residue_names:
            for atom in getattr(molecule, "atoms", ()):
                metadata = getattr(atom, "metadata", None)
                if isinstance(metadata, dict) and _atom_force_field_domain(atom) is None:
                    metadata.setdefault("force_field_domain", "sage")
            molecule_domains.add("sage")
        if len(molecule_domains) == 1 and not hasattr(molecule, "force_field_domain"):
            domain_value = next(iter(molecule_domains))
            try:
                molecule.force_field_domain = domain_value
            except AttributeError:
                properties = getattr(molecule, "properties", None)
                if isinstance(properties, dict):
                    properties.setdefault("force_field_domain", domain_value)


def _sage_residue_names_from_config(config: Any) -> set[str]:
    """Return disconnected Sage residue names declared by solvent and free polymer config."""
    names: set[str] = set()
    solvent = getattr(config, "solvent", None)
    for cosolvent in getattr(solvent, "co_solvents", ()) or ():
        residue_name = getattr(cosolvent, "residue_name", None) or getattr(cosolvent, "name", None)
        if residue_name:
            names.add(str(residue_name).strip().upper()[:3])
    polymers = getattr(config, "polymers", None)
    if polymers is not None and getattr(polymers, "enabled", False):
        for monomer in getattr(polymers, "monomers", ()) or ():
            residue_name = getattr(monomer, "residue_name", None) or getattr(monomer, "name", None)
            if residue_name:
                names.add(str(residue_name).strip().upper()[:3])
    return names


def _register_disconnected_sage_template_generator(
    force_field: Any, topology: Any, config: Any, construction: Any
) -> tuple[dict[str, Any], ...]:
    """Register complete disconnected precharged Sage-domain molecules.

    Only molecules with explicit Sage-domain provenance and assigned partial
    charges are admitted. Strict glycan molecules and the product protein
    component are excluded so GLYCAM/NLN ownership remains authoritative. A Sage
    molecule covalently connected to protein or glycan atoms fails because this
    implementation cannot prove cross-boundary bonded and exception provenance.
    """
    sage_molecules = []
    diagnostics = []
    sage_residue_names = _sage_residue_names_from_config(config)
    for molecule_index, molecule in enumerate(getattr(topology, "molecules", ())):
        domain = _molecule_force_field_domain(molecule)
        if domain is None and _molecule_matches_sage_residue_name(molecule, sage_residue_names):
            domain = "sage"
        if domain == "glycan":
            diagnostics.append(
                {"molecule_index": molecule_index, "domain": domain, "routed_to_sage": False}
            )
            continue
        if domain != "sage":
            continue
        if _sage_molecule_has_amber_glycam_boundary_bond(molecule):
            raise ValueError(
                "Covalently attached Sage-domain molecule crosses an Amber/GLYCAM boundary; "
                "this route is unsupported because cross-boundary bonded and exception "
                "provenance is not audited"
            )
        if not _molecule_has_assigned_partial_charges(molecule):
            raise ValueError(
                f"Sage-domain molecule {molecule_index} is missing trusted assigned partial charges"
            )
        sage_molecules.append(molecule)
        diagnostics.append(
            {
                "molecule_index": molecule_index,
                "component_id": _molecule_component_id(molecule, molecule_index),
                "domain": "sage",
                "routed_to_sage": True,
                "source_residue_segments": _molecule_source_residue_segments(molecule),
            }
        )
    if not sage_molecules:
        return tuple(diagnostics)
    try:
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "Disconnected Sage components require openmmforcefields and "
            "SMIRNOFFTemplateGenerator in the PolyzyMD pixi environment"
        ) from exc

    small_molecule_force_field = getattr(
        getattr(config, "force_field", None), "small_molecule", None
    )
    generator = SMIRNOFFTemplateGenerator(
        molecules=sage_molecules,
        forcefield=small_molecule_force_field or "openff-2.0.0.offxml",
    )
    force_field.registerTemplateGenerator(generator.generator)
    version = _openmmforcefields_version()
    return tuple({**entry, "openmmforcefields_version": version} for entry in diagnostics)


def _molecule_force_field_domain(molecule: Any) -> str | None:
    """Return the force-field domain recorded on a molecule or its atoms."""
    for key in ("force_field_domain", "ff_domain", "domain"):
        value = getattr(molecule, key, None)
        if value is not None:
            return str(value).strip().lower()
    properties = getattr(molecule, "properties", {}) or {}
    if isinstance(properties, dict):
        for key in ("force_field_domain", "ff_domain", "domain"):
            value = properties.get(key)
            if value is not None:
                return str(value).strip().lower()
    atom_domains = {
        str((getattr(atom, "metadata", {}) or {}).get("force_field_domain", "")).strip().lower()
        for atom in getattr(molecule, "atoms", ())
    }
    atom_domains.discard("")
    if len(atom_domains) == 1:
        return atom_domains.pop()
    if len(atom_domains) > 1:
        raise ValueError(f"Molecule has ambiguous force-field domains: {sorted(atom_domains)}")
    return None


def _sage_molecule_has_amber_glycam_boundary_bond(molecule: Any) -> bool:
    """Return whether a Sage molecule has an explicit Amber/GLYCAM boundary bond."""
    if _molecule_is_covalently_connected_to_strict_domain(molecule):
        return True
    strict_domains = {"protein", "amber", "glycan", "glycam", "protein_modified_nln"}
    for bond in getattr(molecule, "bonds", ()) or ():
        domain1 = _atom_force_field_domain(getattr(bond, "atom1", None))
        domain2 = _atom_force_field_domain(getattr(bond, "atom2", None))
        if {domain1, domain2} & strict_domains and "sage" in {domain1, domain2}:
            return True
    return False


def _atom_force_field_domain(atom: Any) -> str | None:
    """Return normalized force-field domain metadata from an atom-like object."""
    metadata = getattr(atom, "metadata", {}) or {}
    value = metadata.get("force_field_domain")
    if value is None:
        residue_name = str(metadata.get("residue_name", "") or "").strip().upper()
        if residue_name == "ASX":
            return "protein_modified_nln"
        if residue_name[:1].isdigit():
            return "glycan"
        return None
    return str(value).strip().lower() or None


def _molecule_matches_sage_residue_name(molecule: Any, sage_residue_names: set[str]) -> bool:
    """Return whether a molecule appears to be a configured disconnected Sage residue."""
    if not sage_residue_names:
        return False
    molecule_names = {
        str(getattr(molecule, key, "") or "").strip().upper()[:3]
        for key in ("name", "component_id")
    }
    if molecule_names & sage_residue_names:
        return True
    for atom in getattr(molecule, "atoms", ()):
        metadata = getattr(atom, "metadata", {}) or {}
        residue_name = str(metadata.get("residue_name", "") or "").strip().upper()[:3]
        if residue_name in sage_residue_names:
            return True
    return False


def _molecule_has_assigned_partial_charges(molecule: Any) -> bool:
    """Return whether a molecule carries complete trusted partial charges."""
    charges = getattr(molecule, "partial_charges", None)
    if charges is None:
        charges = [getattr(atom, "partial_charge", None) for atom in getattr(molecule, "atoms", ())]
    if charges is None:
        return False
    try:
        values = list(charges)
    except TypeError:
        return False
    return bool(values) and all(value is not None for value in values)


def _molecule_is_covalently_connected_to_strict_domain(molecule: Any) -> bool:
    """Return whether a Sage molecule declares cross-domain covalent provenance."""
    metadata = getattr(molecule, "metadata", {}) or {}
    if metadata.get("covalently_connected_to") in {"protein", "glycan", "amber", "glycam"}:
        return True
    return bool(metadata.get("cross_force_field_boundary", False))


def _molecule_component_id(molecule: Any, molecule_index: int) -> str:
    """Return a stable component identifier for Sage-domain audit entries."""
    metadata = getattr(molecule, "metadata", {}) or {}
    return str(metadata.get("component_id") or metadata.get("name") or f"molecule-{molecule_index}")


def _molecule_source_residue_segments(molecule: Any) -> tuple[dict[str, Any], ...]:
    """Return original residue segmentation provenance for a molecule."""
    segments: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for atom_index, atom in enumerate(getattr(molecule, "atoms", ()) or ()):  # pragma: no branch
        identity = _atom_identity(atom)
        key = (
            identity["chain_id"],
            identity["residue_name"],
            identity["residue_number"],
            identity["insertion_code"],
        )
        segment = segments.setdefault(
            key,
            {
                "chain_id": identity["chain_id"],
                "residue_name": identity["residue_name"],
                "residue_number": identity["residue_number"],
                "insertion_code": identity["insertion_code"],
                "atom_indices": [],
                "atom_names": [],
            },
        )
        segment["atom_indices"].append(atom_index)
        segment["atom_names"].append(identity["atom_name"])
    return tuple(segments.values())


def _openmmforcefields_version() -> str | None:
    """Return the installed openmmforcefields version when available."""
    try:
        import openmmforcefields
    except ImportError:  # pragma: no cover - environment dependent
        return None
    return str(getattr(openmmforcefields, "__version__", "unknown"))


def _preflight_glycam_template_matches(
    force_field: Any, topology: Any
) -> tuple[dict[str, Any], ...]:
    """Use OpenMM template matching to verify GLYCAM/NLN residues before system creation."""
    if not hasattr(force_field, "getMatchingTemplates"):
        return ()
    residues = list(topology.residues())
    try:
        matches = force_field.getMatchingTemplates(topology, ignoreExternalBonds=False)
    except TypeError:
        matches = force_field.getMatchingTemplates(topology)
    except Exception as exc:  # noqa: BLE001 - OpenMM reports template failures this way
        blocked_residue_name = _template_error_residue_name(exc)
        if blocked_residue_name is not None and not _strict_glycam_residue_name(
            blocked_residue_name
        ):
            return _strict_glycam_preflight_blocked_by_non_strict_residue(residues, exc)
        raise ValueError(f"Strict GLYCAM preflight failed during template matching: {exc}") from exc
    diagnostics = []
    for residue, match in zip(residues, matches, strict=False):
        template_name = getattr(match, "name", str(match)) if match is not None else None
        if match is None and _is_strict_glycam_residue(residue):
            raise ValueError(
                f"Strict GLYCAM preflight could not match residue {_residue_label(residue)}; "
                "check residue names, atom names, link atoms, and GLYCAM template coverage"
            )
        if _is_strict_glycam_residue(residue):
            diagnostics.append(
                {
                    "residue": _residue_label(residue),
                    "template": template_name,
                    "domain": "glycan" if residue.name != "ASX" else "protein_modified_nln",
                }
            )
    return tuple(diagnostics)


def _template_error_residue_name(exc: Exception) -> str | None:
    """Return the residue name reported by an OpenMM template-match error."""
    match = re.search(r"resname <([^>]+)>", str(exc))
    if match is None:
        return None
    return match.group(1).strip().upper()


def _strict_glycam_preflight_blocked_by_non_strict_residue(
    residues: list[Any], exc: Exception
) -> tuple[dict[str, Any], ...]:
    """Return strict-residue diagnostics when Sage residues block global matching."""
    reason = str(exc).splitlines()[0]
    return tuple(
        {
            "residue": _residue_label(residue),
            "template": "not_evaluated_non_strict_residue_blocked_global_preflight",
            "domain": "glycan" if residue.name != "ASX" else "protein_modified_nln",
            "blocked_by": reason,
        }
        for residue in residues
        if _is_strict_glycam_residue(residue)
    )


def _is_strict_glycam_residue(residue: Any) -> bool:
    """Return whether a residue is owned by strict GLYCAM/NLN routing."""
    if _strict_glycam_residue_name(residue.name):
        return True
    domain = getattr(residue, "force_field_domain", None)
    if domain is not None:
        return str(domain).strip().lower() == "glycan"
    return (
        residue.name not in SUPPORTED_WATER_RESIDUES | SUPPORTED_ION_RESIDUES
        and residue.name[:1].isdigit()
    )


def _strict_glycam_residue_name(residue_name: str) -> bool:
    """Return whether a residue name is intrinsically strict GLYCAM/NLN."""
    residue_name = residue_name.strip().upper()
    return residue_name == "ASX" or residue_name[:1].isdigit()


def _openff_topology_to_openmm_for_glycam(
    topology: Any,
    construction: Any,
    *,
    sage_residue_names: set[str] | None = None,
) -> _ConvertedTopology:
    """Convert an authoritative OpenFF topology graph into an OpenMM topology."""
    from openmm import Vec3, unit
    from openmm.app import Element, Topology

    positions_nm = _positions_nm(topology)
    omm_topology = Topology()
    atom_map: dict[Any, Any] = {}
    residue_map: dict[tuple[str, str, str, str], Any] = {}
    renamed_atoms: list[dict[str, Any]] = []
    sage_template_units: list[dict[str, Any]] = []

    for molecule_index, molecule in enumerate(topology.molecules):
        chain_by_id: dict[str, Any] = {}
        molecule_domain = _molecule_force_field_domain(molecule)
        if molecule_domain is None and _molecule_matches_sage_residue_name(
            molecule, sage_residue_names or set()
        ):
            molecule_domain = "sage"
        if molecule_domain == "sage":
            if _sage_molecule_has_amber_glycam_boundary_bond(molecule):
                raise ValueError(
                    "Covalently attached Sage-domain molecule crosses an Amber/GLYCAM boundary; "
                    "this route is unsupported because cross-boundary bonded and exception "
                    "provenance is not audited"
                )
            template_unit = _add_disconnected_sage_template_unit(
                omm_topology,
                molecule,
                molecule_index=molecule_index,
                atom_map=atom_map,
            )
            sage_template_units.append(template_unit)
            continue
        for atom in molecule.atoms:
            identity = _atom_identity(atom)
            if _is_asx_hd22(identity):
                renamed_atoms.append({**identity, "native_atom_name": "HD21"})
                identity = {**identity, "atom_name": "HD21"}
            chain = chain_by_id.get(identity["chain_id"])
            if chain is None:
                chain = omm_topology.addChain(identity["chain_id"])
                chain_by_id[identity["chain_id"]] = chain
            residue_key = (
                str(molecule_index),
                identity["chain_id"],
                identity["residue_name"],
                identity["residue_number"],
                identity["insertion_code"],
            )
            residue = residue_map.get(residue_key)
            if residue is None:
                residue = omm_topology.addResidue(
                    identity["residue_name"],
                    chain,
                    id=identity["residue_number"],
                    insertionCode=identity["insertion_code"],
                )
                _copy_residue_domain_metadata(residue, atom)
                residue_map[residue_key] = residue
            element = Element.getByAtomicNumber(int(atom.atomic_number))
            atom_map[atom] = omm_topology.addAtom(identity["atom_name"], element, residue)

    kept_positions = []
    topology_atom_index = 0
    for molecule in topology.molecules:
        for atom in molecule.atoms:
            if atom not in atom_map:
                topology_atom_index += 1
                continue
            kept_positions.append(positions_nm[topology_atom_index])
            topology_atom_index += 1

    for molecule in topology.molecules:
        for bond in molecule.bonds:
            atom1 = atom_map.get(bond.atom1)
            atom2 = atom_map.get(bond.atom2)
            if atom1 is None or atom2 is None:
                continue
            omm_topology.addBond(atom1, atom2)
    _set_box_vectors(omm_topology, topology)

    asx_residues = _require_asx_residues(omm_topology)
    _require_asx_canonical_atoms(asx_residues, renamed_atoms)
    crosslink_pairs = _require_crosslinks(omm_topology, construction)
    positions = [Vec3(float(x), float(y), float(z)) for x, y, z in kept_positions] * unit.nanometer
    return _ConvertedTopology(
        topology=omm_topology,
        positions=positions,
        asx_residues=asx_residues,
        crosslink_pairs=crosslink_pairs,
        renamed_atoms=tuple(renamed_atoms),
        sage_template_units=tuple(sage_template_units),
    )


def _add_disconnected_sage_template_unit(
    omm_topology: Any,
    molecule: Any,
    *,
    molecule_index: int,
    atom_map: dict[Any, Any],
) -> dict[str, Any]:
    """Add one OpenMM residue containing a complete disconnected Sage molecule.

    This one-residue representation is only a SMIRNOFF template-matching
    implementation detail. The source residue segmentation is preserved in audit
    provenance so polymer monomer identities remain recoverable.
    """
    from openmm.app import Element

    atoms = tuple(getattr(molecule, "atoms", ()) or ())
    if not atoms:
        raise ValueError(f"Sage-domain molecule {molecule_index} contains no atoms")
    first_identity = _atom_identity(atoms[0])
    chain = omm_topology.addChain(first_identity["chain_id"])
    residue_name = _sage_template_residue_name(molecule, first_identity)
    residue = omm_topology.addResidue(
        residue_name,
        chain,
        id=_sage_template_residue_id(molecule_index, first_identity),
        insertionCode="",
    )
    try:
        residue.force_field_domain = "sage"
        residue.source_residue_segments = _molecule_source_residue_segments(molecule)
        residue.source_component_id = _molecule_component_id(molecule, molecule_index)
    except AttributeError:
        pass
    for atom in atoms:
        identity = _atom_identity(atom)
        element = Element.getByAtomicNumber(int(atom.atomic_number))
        atom_map[atom] = omm_topology.addAtom(identity["atom_name"], element, residue)
    return {
        "molecule_index": molecule_index,
        "component_id": _molecule_component_id(molecule, molecule_index),
        "template_residue": _residue_label(residue),
        "template_residue_name": residue_name,
        "template_residue_id": residue.id,
        "atom_count": len(atoms),
        "source_residue_segments": _molecule_source_residue_segments(molecule),
    }


def _sage_template_residue_name(molecule: Any, first_identity: dict[str, Any]) -> str:
    """Return a deterministic residue name for a complete Sage template unit."""
    metadata = getattr(molecule, "metadata", {}) or {}
    for value in (
        metadata.get("template_residue_name"),
        metadata.get("residue_name"),
        first_identity["residue_name"],
    ):
        normalized = str(value or "").strip().upper()
        if normalized:
            return normalized[:3]
    return "SGE"


def _sage_template_residue_id(molecule_index: int, first_identity: dict[str, Any]) -> str:
    """Return a deterministic unique residue id for a collapsed Sage molecule."""
    _ = first_identity
    return f"S{molecule_index + 1}"


def _positions_nm(topology: Any) -> np.ndarray:
    """Return OpenFF topology positions as a nanometer array."""
    positions = topology.get_positions()
    if hasattr(positions, "m_as"):
        return np.asarray(positions.m_as("nanometer"), dtype=float)
    return np.asarray(positions, dtype=float)


def _atom_identity(atom: Any) -> dict[str, Any]:
    """Return canonical PDB identity metadata for an OpenFF atom."""
    metadata = getattr(atom, "metadata", {}) or {}
    residue_number = str(metadata.get("residue_number", metadata.get("residue_id", "1"))).strip()
    return {
        "chain_id": str(metadata.get("chain_id", "A")).strip() or "A",
        "residue_name": str(metadata.get("residue_name", metadata.get("residue", "UNK")))
        .strip()
        .upper(),
        "residue_number": residue_number or "1",
        "insertion_code": str(metadata.get("insertion_code", "") or "").strip(),
        "atom_name": str(metadata.get("atom_name", getattr(atom, "name", ""))).strip().upper(),
        "atomic_number": int(atom.atomic_number),
    }


def _copy_residue_domain_metadata(residue: Any, source_atom: Any) -> None:
    """Copy force-field domain provenance from OpenFF atom metadata to an OpenMM residue."""
    metadata = getattr(source_atom, "metadata", {}) or {}
    domain = metadata.get("force_field_domain")
    if domain is None:
        return
    try:
        residue.force_field_domain = str(domain).strip().lower()
    except AttributeError:
        pass


def _is_asx_hd22(identity: dict[str, Any]) -> bool:
    """Return whether the identity is the ASX HD22 atom renamed for NLN."""
    return identity["residue_name"] == "ASX" and identity["atom_name"] == "HD22"


def _set_box_vectors(omm_topology: Any, topology: Any) -> None:
    """Copy OpenFF box vectors onto the OpenMM topology."""
    box_vectors = getattr(topology, "box_vectors", None)
    if box_vectors is None:
        return
    from openmm import Vec3, unit

    if hasattr(box_vectors, "m_as"):
        matrix = np.asarray(box_vectors.m_as("nanometer"), dtype=float)
    else:
        matrix = np.asarray(box_vectors, dtype=float)
    omm_topology.setPeriodicBoxVectors(tuple(Vec3(*row) for row in matrix) * unit.nanometer)


def _require_asx_residues(omm_topology: Any) -> tuple[Any, ...]:
    """Return ASX residues that will be mapped to GLYCAM NLN templates."""
    residues = [residue for residue in omm_topology.residues() if residue.name == "ASX"]
    if not residues:
        raise ValueError("Strict native GLYCAM mode requires at least one modified ASX residue")
    return tuple(residues)


def _require_asx_canonical_atoms(
    asx_residues: tuple[Any, ...], renamed_atoms: list[dict[str, Any]]
) -> None:
    """Validate ASX atom names required for native GLYCAM NLN mapping."""
    expected_renames = 0
    for asx_residue in asx_residues:
        atom_name_list = [atom.name for atom in asx_residue.atoms()]
        atom_names = set(atom_name_list)
        if "ND2" not in atom_names:
            raise ValueError(
                f"Native GLYCAM ASX residue {_residue_label(asx_residue)} must use atom ND2"
            )
        retained_hydrogens = atom_names.intersection({"HD21", "HD22"})
        expected_renames += 1
        if retained_hydrogens != {"HD21"} or atom_name_list.count("HD21") != 1:
            raise ValueError(
                "Native GLYCAM ASX residue must contain exactly one retained amide hydrogen: "
                f"source HD22 renamed to native HD21 at {_residue_label(asx_residue)}"
            )
    if len(renamed_atoms) != expected_renames:
        raise ValueError(
            "Native GLYCAM ASX residue must contain exactly one retained amide hydrogen: "
            "source HD22 renamed to native HD21"
        )


def _require_crosslinks(omm_topology: Any, construction: Any) -> tuple[tuple[Any, Any], ...]:
    """Return and validate all ASX ND2 to GLYCAM residue crosslink atoms."""
    _ = construction
    crosslinks = []
    for atom1, atom2 in omm_topology.bonds():
        if _is_asx_link_atom(atom1) and atom2.residue.name != "ASX":
            crosslinks.append((atom1, atom2))
        elif _is_asx_link_atom(atom2) and atom1.residue.name != "ASX":
            crosslinks.append((atom2, atom1))
    if not crosslinks:
        raise ValueError("Native GLYCAM topology is missing explicit ASX ND2--glycan link bond")
    return tuple(crosslinks)


def _is_asx_link_atom(atom: Any) -> bool:
    """Return whether an atom is the modified ASN/NLN linkage nitrogen."""
    return atom.residue.name == "ASX" and atom.name == "ND2"


def _single_atom(atoms: list[Any], *, residue_name: str, atom_name: str) -> Any:
    """Return one atom selected by residue and atom names."""
    matches = [
        atom for atom in atoms if atom.residue.name == residue_name and atom.name == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Native GLYCAM mode requires exactly one {residue_name} {atom_name}; "
            f"found {len(matches)}"
        )
    return matches[0]


def _build_native_glycam_audit(
    topology: Any,
    system: Any,
    asx_residues: tuple[Any, ...],
    crosslink_pairs: tuple[tuple[Any, Any], ...],
    *,
    renamed_atoms: tuple[dict[str, Any], ...],
    template_matches: tuple[dict[str, Any], ...] = (),
    sage_components: tuple[dict[str, Any], ...] = (),
    construction: Any = None,
) -> dict[str, Any]:
    """Build JSON-safe evidence for native GLYCAM parameterization."""
    import openmm
    from openmm import unit

    atoms = list(topology.atoms())
    residues = list(topology.residues())
    residue_counts = Counter(residue.name for residue in residues)
    nonbonded = _nonbonded_force(system)
    linkage_audits = []
    for linkage_index, (nd2, glycan_atom) in enumerate(crosslink_pairs):
        bonded_terms = _local_bonded_terms(system, atoms, nd2.index, glycan_atom.index)
        exceptions = _local_nonbonded_exceptions(
            topology, nonbonded, atoms, nd2.index, glycan_atom.index
        )
        linkage_audits.append(
            {
                "index": linkage_index,
                "identity": f"{nd2.residue.name} {nd2.name}--{glycan_atom.residue.name} {glycan_atom.name}",
                "atoms": [_atom_label(nd2), _atom_label(glycan_atom)],
                "bond_present": True,
                "adjacent_terms": bonded_terms,
                "local_exclusions_and_14_exceptions": exceptions,
            }
        )
    charges = [
        float(nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge))
        for index in range(nonbonded.getNumParticles())
    ]
    return {
        "route": "native_openmm_glycam",
        "route_invariants": dict(NATIVE_GLYCAM_ROUTE_INVARIANTS),
        "openmm_version": openmm.version.version,
        "xml_stack": list(GLYCAM_XML_STACK),
        "domain_assignments": _domain_assignment_audit(residues, sage_components),
        "glycam_template_matches": tuple(template_matches),
        "sage_template_generator": {
            "components": tuple(entry for entry in sage_components if entry.get("routed_to_sage")),
            "proof_no_glycan_entered_sage": all(
                not entry.get("routed_to_sage")
                for entry in sage_components
                if entry.get("domain") == "glycan"
            ),
        },
        "unsupported_boundary_diagnostics": {
            "covalently_attached_sage_polymer_across_amber_glycam_boundary": "unsupported",
            "reason": (
                "No audited force-field provenance defines cross-boundary bonded terms and "
                "NonbondedForce exceptions for covalently attached Sage plus GLYCAM domains."
            ),
        },
        "residue_templates": {_residue_label(asx_residue): "NLN" for asx_residue in asx_residues},
        "counts": {
            "atoms": topology.getNumAtoms(),
            "residues": len(residues),
            "bonds": topology.getNumBonds(),
            "particles": system.getNumParticles(),
            "constraints": system.getNumConstraints(),
            "residue_names": dict(sorted(residue_counts.items())),
        },
        "total_charge_e": sum(charges),
        "pme_settings": _pme_settings(nonbonded),
        "crosslinks": tuple(
            {
                "identity": linkage["identity"],
                "atoms": linkage["atoms"],
                "bond_present": linkage["bond_present"],
            }
            for linkage in linkage_audits
        ),
        "attachment_provenance": _attachment_provenance_audit(construction),
        "nln_hydrogen_policy": {
            "source_retained": "HD22",
            "native_retained": "HD21",
            "renamed": tuple(renamed_atoms),
        },
        "linkages": tuple(linkage_audits),
        "adjacent_terms": linkage_audits[0]["adjacent_terms"] if linkage_audits else {},
        "local_exclusions_and_14_exceptions": (
            linkage_audits[0]["local_exclusions_and_14_exceptions"] if linkage_audits else ()
        ),
    }


def _domain_assignment_audit(
    residues: list[Any], sage_components: tuple[dict[str, Any], ...]
) -> dict[str, Any]:
    """Return residue and component domain provenance for exact export."""
    return {
        "residues": tuple(
            {
                "residue": _residue_label(residue),
                "domain": _residue_domain(residue),
                "source_force_field": _residue_source_force_field(residue),
                **_residue_source_provenance(residue),
            }
            for residue in residues
        ),
        "sage_components": sage_components,
    }


def _residue_domain(residue: Any) -> str:
    """Return the force-field domain label for an OpenMM residue."""
    domain = getattr(residue, "force_field_domain", None)
    if domain is not None and str(domain).strip().lower() == "sage":
        return "sage"
    if residue.name == "ASX":
        return "protein_modified_nln"
    if residue.name in SUPPORTED_WATER_RESIDUES:
        return "water"
    if residue.name in SUPPORTED_ION_RESIDUES:
        return "ion"
    if _is_strict_glycam_residue(residue):
        return "glycan"
    return "protein"


def _residue_source_provenance(residue: Any) -> dict[str, Any]:
    """Return optional source provenance attached during topology conversion."""
    provenance: dict[str, Any] = {}
    source_segments = getattr(residue, "source_residue_segments", None)
    if source_segments is not None:
        provenance["source_residue_segments"] = source_segments
    source_component_id = getattr(residue, "source_component_id", None)
    if source_component_id is not None:
        provenance["source_component_id"] = source_component_id
    return provenance


def _residue_source_force_field(residue: Any) -> str:
    """Return the source force field for an OpenMM residue domain."""
    domain = _residue_domain(residue)
    if domain == "protein":
        return "Amber ff14SB"
    if domain == "protein_modified_nln":
        return "GLYCAM06 NLN"
    if domain == "glycan":
        return "GLYCAM06"
    if domain == "water":
        return "Amber TIP3P"
    if domain == "ion":
        return "Amber/TIP3P standard ion template"
    return "OpenFF Sage"


def _attachment_provenance_audit(construction: Any) -> tuple[dict[str, Any], ...]:
    """Return lightweight attachment provenance entries when available."""
    attachments = getattr(construction, "attachment_specs", None) or getattr(
        construction, "attachments", ()
    )
    entries = []
    for index, attachment in enumerate(attachments or ()):  # pragma: no branch - tiny audit loop
        entries.append(
            {
                "index": index,
                "name": str(getattr(attachment, "name", f"attachment-{index}")),
                "force_field_domain": str(
                    getattr(getattr(attachment, "moiety", None), "force_field_domain", "glycan")
                ),
            }
        )
    return tuple(entries)


def _nonbonded_force(system: Any) -> Any:
    """Return the system NonbondedForce or fail clearly."""
    for force in system.getForces():
        if force.__class__.__name__ == "NonbondedForce":
            return force
    raise RuntimeError("Native GLYCAM audit could not find OpenMM NonbondedForce")


def _pme_settings(nonbonded: Any) -> dict[str, Any]:
    """Return JSON-safe PME settings from a NonbondedForce."""
    from openmm import unit

    return {
        "method": int(nonbonded.getNonbondedMethod()),
        "cutoff_nm": float(nonbonded.getCutoffDistance().value_in_unit(unit.nanometer)),
        "ewald_error_tolerance": float(nonbonded.getEwaldErrorTolerance()),
        "dispersion_correction": bool(nonbonded.getUseDispersionCorrection()),
    }


def _local_bonded_terms(
    system: Any, atoms: list[Any], nd2_index: int, c1_index: int
) -> dict[str, Any]:
    """Collect bond, angle, and torsion terms adjacent to the crosslink."""
    from openmm import unit

    terms: dict[str, list[dict[str, Any]]] = {"bonds": [], "angles": [], "torsions": []}
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "HarmonicBondForce":
            for index in range(force.getNumBonds()):
                i, j, length, k_value = force.getBondParameters(index)
                if {i, j} == {nd2_index, c1_index}:
                    terms["bonds"].append(
                        {
                            "index": index,
                            "category": "crosslink_bond",
                            "atoms": [_atom_label(atoms[i]), _atom_label(atoms[j])],
                            "length_nm": float(length.value_in_unit(unit.nanometer)),
                            "k_kj_mol_nm2": float(
                                k_value.value_in_unit(unit.kilojoule_per_mole / unit.nanometer**2)
                            ),
                        }
                    )
        elif name == "HarmonicAngleForce":
            for index in range(force.getNumAngles()):
                i, j, k_index, angle, k_value = force.getAngleParameters(index)
                if nd2_index in {i, j, k_index} and c1_index in {i, j, k_index}:
                    category = _angle_category(i, j, k_index, nd2_index, c1_index)
                    terms["angles"].append(
                        {
                            "index": index,
                            "category": category,
                            "atoms": [_atom_label(atoms[x]) for x in (i, j, k_index)],
                            "angle_rad": float(angle.value_in_unit(unit.radian)),
                            "k_kj_mol_rad2": float(
                                k_value.value_in_unit(unit.kilojoule_per_mole / unit.radian**2)
                            ),
                        }
                    )
        elif name == "PeriodicTorsionForce":
            for index in range(force.getNumTorsions()):
                i, j, k_index, l_index, periodicity, phase, k_value = force.getTorsionParameters(
                    index
                )
                if nd2_index in {i, j, k_index, l_index} and c1_index in {
                    i,
                    j,
                    k_index,
                    l_index,
                }:
                    category = _torsion_category(i, j, k_index, l_index, nd2_index, c1_index)
                    terms["torsions"].append(
                        {
                            "index": index,
                            "category": category,
                            "atoms": [_atom_label(atoms[x]) for x in (i, j, k_index, l_index)],
                            "periodicity": int(periodicity),
                            "phase_rad": float(phase.value_in_unit(unit.radian)),
                            "k_kj_mol": float(k_value.value_in_unit(unit.kilojoule_per_mole)),
                        }
                    )
    return terms


def _local_nonbonded_exceptions(
    topology: Any, nonbonded: Any, atoms: list[Any], nd2_index: int, c1_index: int
) -> list[dict[str, Any]]:
    """Collect exclusions and 1-4 exceptions within four bonds of the linkage."""
    from openmm import unit

    adjacency = _topology_adjacency(topology)
    graph_distance = _graph_distances(adjacency, (nd2_index, c1_index), max_distance=4)
    exceptions = []
    for index in range(nonbonded.getNumExceptions()):
        i, j, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(index)
        if i not in graph_distance or j not in graph_distance:
            continue
        pair_distance = _shortest_pair_distance(adjacency, i, j, max_distance=4)
        chargeprod_e2 = float(chargeprod.value_in_unit(unit.elementary_charge**2))
        epsilon_kj_mol = float(epsilon.value_in_unit(unit.kilojoule_per_mole))
        exceptions.append(
            {
                "index": index,
                "category": _exception_category(pair_distance, chargeprod_e2, epsilon_kj_mol),
                "atoms": [_atom_label(atoms[i]), _atom_label(atoms[j])],
                "graph_distance_pair": [graph_distance[i], graph_distance[j]],
                "pair_graph_distance": pair_distance,
                "chargeprod_e2": chargeprod_e2,
                "sigma_nm": float(sigma.value_in_unit(unit.nanometer)),
                "epsilon_kj_mol": epsilon_kj_mol,
            }
        )
    return exceptions


def _topology_adjacency(topology: Any) -> dict[int, list[int]]:
    """Return an atom-index adjacency map for an OpenMM topology."""
    adjacency: dict[int, list[int]] = {}
    for atom1, atom2 in topology.bonds():
        adjacency.setdefault(atom1.index, []).append(atom2.index)
        adjacency.setdefault(atom2.index, []).append(atom1.index)
    return adjacency


def _graph_distances(
    adjacency: dict[int, list[int]], starts: tuple[int, ...], *, max_distance: int
) -> dict[int, int]:
    """Return graph distances from selected atom indices."""
    distances = dict.fromkeys(starts, 0)
    queue = deque(starts)
    while queue:
        index = queue.popleft()
        if distances[index] >= max_distance:
            continue
        for neighbor in adjacency.get(index, []):
            if neighbor in distances:
                continue
            distances[neighbor] = distances[index] + 1
            queue.append(neighbor)
    return distances


def _shortest_pair_distance(
    adjacency: dict[int, list[int]], start: int, target: int, *, max_distance: int
) -> int | None:
    """Return the shortest graph distance between two atom indices."""
    distances = _graph_distances(adjacency, (start,), max_distance=max_distance)
    return distances.get(target)


def _angle_category(i: int, j: int, k_index: int, nd2_index: int, c1_index: int) -> str:
    """Classify a crosslink angle by its central atom."""
    if j == nd2_index and c1_index in {i, k_index}:
        return "asx_side_angle"
    if j == c1_index and nd2_index in {i, k_index}:
        return "glycan_side_angle"
    return "other_crosslink_angle"


def _torsion_category(
    i: int, j: int, k_index: int, l_index: int, nd2_index: int, c1_index: int
) -> str:
    """Classify a crosslink torsion by whether the linkage is the central bond."""
    if (j, k_index) in {(nd2_index, c1_index), (c1_index, nd2_index)}:
        return "proper_crosslink_torsion"
    if {nd2_index, c1_index}.issubset({i, j, k_index, l_index}):
        return "other_crosslink_torsion"
    return "unrelated_torsion"


def _exception_category(
    pair_distance: int | None, chargeprod_e2: float, epsilon_kj_mol: float
) -> str:
    """Classify a local nonbonded exception by graph distance and scaling."""
    is_zero = abs(chargeprod_e2) < 1e-12 and abs(epsilon_kj_mol) < 1e-12
    if pair_distance == 1 and is_zero:
        return "zero_12_exclusion"
    if pair_distance == 2 and is_zero:
        return "zero_13_exclusion"
    if pair_distance == 3 and not is_zero:
        return "scaled_14_exception"
    return "other_local_exception"


def _require_essential_linkage_terms(audit: dict[str, Any]) -> None:
    """Fail if OpenMM did not instantiate essential crosslink terms."""
    linkages = audit.get("linkages") or (
        {
            "index": 0,
            "adjacent_terms": audit["adjacent_terms"],
            "local_exclusions_and_14_exceptions": audit["local_exclusions_and_14_exceptions"],
        },
    )
    for linkage in linkages:
        _require_single_linkage_terms(linkage)


def _require_single_linkage_terms(linkage: dict[str, Any]) -> None:
    """Fail if one linkage lacks essential bond, angle, torsion, or exception terms."""
    prefix = f"linkage {linkage.get('index', 0)} "
    terms = linkage["adjacent_terms"]
    if not terms["bonds"]:
        raise RuntimeError(f"Native GLYCAM audit missing category: {prefix}crosslink_bond")
    angle_categories = {term["category"] for term in terms["angles"]}
    for category in ("asx_side_angle", "glycan_side_angle"):
        if category not in angle_categories:
            raise RuntimeError(f"Native GLYCAM audit missing category: {prefix}{category}")
    torsion_categories = {term["category"] for term in terms["torsions"]}
    if "proper_crosslink_torsion" not in torsion_categories:
        raise RuntimeError(
            f"Native GLYCAM audit missing category: {prefix}proper_crosslink_torsion"
        )
    exception_categories = {
        term["category"] for term in linkage["local_exclusions_and_14_exceptions"]
    }
    for category in ("zero_12_exclusion", "zero_13_exclusion", "scaled_14_exception"):
        if category not in exception_categories:
            raise RuntimeError(f"Native GLYCAM audit missing category: {prefix}{category}")


def _atom_label(atom: Any) -> str:
    """Return a stable human-readable atom label."""
    return f"{_residue_label(atom.residue)}:{atom.name}#{atom.index}"


def _residue_label(residue: Any) -> str:
    """Return a stable human-readable residue label."""
    insertion = getattr(residue, "insertionCode", "") or ""
    return f"{residue.chain.id}:{residue.name}{residue.id}{insertion}"
