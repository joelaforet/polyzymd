"""Production charge bridge for final product-state conjugate templates."""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation._linkage import parse_pdb_atom_records
from polyzymd.builders.conjugation.pablo.charge_records import (
    AtomPartialChargeRecord,
    ChargeBridgeReport,
    ResiduePartialChargeRecord,
    validate_unique_atom_records,
)
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    load_combined_smirnoff_force_field,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_PATCH_NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
PATCH_NAGL_MODEL_ENV = "POLYZYMD_CONJUGATE_PATCH_NAGL_MODEL"
_BRIDGE_SOURCE = "production:product-state-local-nagl-charge-bridge"
_PATCH_MAPPED_SMILES = "[CH3:12][NH:13][C:14](=[O:15])[CH3:16]"
_PATCH_NZ_MAP = 13
_PATCH_CARBONYL_C_MAP = 14
_PATCH_CARBONYL_O_MAP = 15


def build_product_state_charge_bridge(
    *,
    product_state_pablo_library: Any,
    product_topology: Any,
    product_pdb: Path | str,
    source_protein_pdb: Path | str,
    specs: Sequence[Any],
    output_dir: Path | str,
    settings: InterchangeParameterizationSettings | None = None,
    total_charge_tolerance: float = 1.0e-4,
) -> Any:
    """Attach production partial-charge records to a product-state Pablo library.

    The bridge keeps the OpenFF charge template at molecule scope while assigning
    atom charges from local source classes: ff14SB protein atoms, charged polymer
    source templates, and a local NAGL product-state linkage patch.

    Parameters
    ----------
    product_state_pablo_library : Any
        Product-state Pablo library to copy and augment.
    product_topology : Any
        Pablo/OpenFF topology for the crosslinked product before solvation.
    product_pdb : pathlib.Path or str
        Crosslinked product PDB used for residue identity mapping.
    source_protein_pdb : pathlib.Path or str
        Prepared unmodified protein PDB used for canonical ff14SB charges.
    specs : sequence of Any
        Attachment specs with resolved plans and generated fragments.
    output_dir : pathlib.Path or str
        Conjugate construction artifact directory.
    settings : InterchangeParameterizationSettings or None, optional
        Force-field settings for ff14SB source extraction, by default ``None``.
    total_charge_tolerance : float, optional
        Allowed total-charge mismatch before normalization, by default ``1e-4``.

    Returns
    -------
    Any
        A copied product-state Pablo library with ``residue_partial_charges``,
        ``charge_bridge_report``, and ``charge_bridge_report_path`` fields set.
    """
    if total_charge_tolerance < 0:
        raise ValueError("total_charge_tolerance must be non-negative")
    if not specs:
        raise ValueError("Product-state charge bridge requires at least one attachment spec")

    artifact_dir = Path(output_dir)
    product_atoms = tuple(parse_pdb_atom_records(Path(product_pdb)))
    product_names = _product_residue_names(product_state_pablo_library)
    target_molecule = _product_conjugate_molecule(
        product_topology,
        product_names=product_names,
        product_atoms=product_atoms,
    )
    target_identities = _target_identities(target_molecule, product_atoms=product_atoms)
    formal_total = sum(
        _formal_charge_value(getattr(atom, "formal_charge", 0)) for atom in target_molecule.atoms
    )

    records: dict[tuple[str, str, int | None, str, str], AtomPartialChargeRecord] = {}
    diagnostics: list[str] = []

    protein_records = _protein_ff14sb_records(
        product_atoms=product_atoms,
        target_identities=target_identities,
        source_protein_pdb=Path(source_protein_pdb),
        settings=settings,
    )
    _merge_records(records, protein_records)

    polymer_records = _polymer_template_records(specs)
    _merge_records(records, polymer_records)

    patch_records, nagl_model = _local_nagl_patch_records(specs)
    _merge_records(records, patch_records, replace=True)

    missing = [identity for identity in target_identities if identity not in records]
    if missing:
        preview = "; ".join(_format_identity(identity) for identity in missing[:12])
        suffix = "" if len(missing) <= 12 else f"; ... {len(missing) - 12} more"
        raise ValueError(
            "Product-state charge bridge could not assign production charges for final "
            f"conjugate atom(s): {preview}{suffix}"
        )

    correction = formal_total - sum(record.charge_e for record in records.values())
    if abs(correction) > total_charge_tolerance:
        correction_key = _normalization_target(records)
        if correction_key is None:
            raise ValueError(
                "Product-state charge bridge has no polymer/template atom available for "
                f"normalizing total charge mismatch of {correction:.8f} e"
            )
        record = records[correction_key]
        records[correction_key] = record.model_copy(
            update={
                "charge_e": record.charge_e + correction,
                "source": f"{record.source}; normalized by {correction:.8f} e",
                "source_role": "normalization",
            }
        )
        diagnostics.append(
            "Applied a uniform-model total-charge closure correction to one attached-polymer "
            f"template atom: {correction:.8f} e"
        )
    else:
        correction = 0.0

    atom_records = tuple(records[identity] for identity in target_identities)
    validate_unique_atom_records(atom_records)
    residue_records = ResiduePartialChargeRecord.from_atom_records(atom_records)
    report = ChargeBridgeReport(
        success=True,
        source=_BRIDGE_SOURCE,
        nagl_model=nagl_model,
        ff14sb_atom_count=sum(record.source_role == "protein_ff14sb" for record in atom_records),
        polymer_template_atom_count=sum(
            record.source_role == "polymer_template" for record in atom_records
        ),
        local_nagl_patch_atom_count=sum(
            record.source_role == "local_nagl_patch" for record in atom_records
        ),
        normalization_correction_e=correction,
        total_charge_e=sum(record.charge_e for record in atom_records),
        formal_charge_e=formal_total,
        diagnostics=tuple(diagnostics),
        assumptions=(
            "Canonical protein atoms retain ff14SB-style charges from the prepared source protein.",
            "Attached polymer atoms retain existing charged polymer/template charges when mapping is stable.",
            "NHS-Lys linkage atoms are overridden by a local product-state NAGL/AshGC patch.",
            "The complete covalent conjugate is emitted as one charged OpenFF Molecule template.",
            "This bridge is not whole-conjugate AM1-BCC and does not use Gasteiger or formal smoke charges.",
        ),
    )
    report_path = report.write_json(artifact_dir / "product_state_charge_bridge.json")
    report = report.model_copy(update={"json_path": report_path})

    LOGGER.info(
        "Product-state charge bridge: %d template(s), %d ff14SB atom(s), %d polymer-template "
        "atom(s), %d local NAGL patch atom(s), correction %.6g e, model %s",
        1,
        report.ff14sb_atom_count,
        report.polymer_template_atom_count,
        report.local_nagl_patch_atom_count,
        report.normalization_correction_e,
        report.nagl_model or "none",
    )

    updates = {
        "residue_partial_charges": residue_records,
        "charge_bridge_report": report,
        "charge_bridge_report_path": report_path,
    }
    if hasattr(product_state_pablo_library, "model_copy"):
        return product_state_pablo_library.model_copy(update=updates)
    for name, value in updates.items():
        setattr(product_state_pablo_library, name, value)
    return product_state_pablo_library


def _product_residue_names(product_state_pablo_library: Any) -> set[str]:
    """Return product-state residue names from a Pablo library."""
    product_names = {
        str(name).strip().upper()
        for name in tuple(getattr(product_state_pablo_library, "residue_names", ()) or ())
        if str(name).strip()
    }
    product_names.update(
        str(getattr(definition, "residue_name", "")).strip().upper()
        for definition in tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
        if str(getattr(definition, "residue_name", "")).strip()
    )
    if not product_names:
        raise ValueError("Product-state charge bridge requires product residue names")
    return product_names


def _product_conjugate_molecule(
    product_topology: Any,
    *,
    product_names: set[str],
    product_atoms: tuple[Any, ...],
) -> Any:
    """Return the topology molecule containing product-state residues."""
    for molecule in tuple(getattr(product_topology, "molecules", ()) or ()):
        if any(_atom_identity(atom)[1] in product_names for atom in getattr(molecule, "atoms", ())):
            return molecule
    molecules = tuple(getattr(product_topology, "molecules", ()) or ())
    if len(molecules) == 1 and any(
        atom.residue_name.upper() in product_names for atom in product_atoms
    ):
        return molecules[0]
    raise ValueError("Could not locate the product-state conjugate molecule in the Pablo topology")


def _target_identities(
    target_molecule: Any,
    *,
    product_atoms: tuple[Any, ...],
) -> tuple[tuple[str, str, int | None, str, str], ...]:
    """Return target atom identities, falling back to product PDB order."""
    atoms = tuple(getattr(target_molecule, "atoms", ()) or ())
    identities = tuple(_atom_identity(atom) for atom in atoms)
    if all(identity[1] and identity[4] for identity in identities):
        return identities
    if len(atoms) == len(product_atoms):
        return tuple(
            (
                atom.chain_id.strip(),
                atom.residue_name.strip().upper(),
                atom.residue_number,
                atom.insertion_code.strip(),
                atom.atom_name.strip(),
            )
            for atom in product_atoms
        )
    return identities


def _protein_ff14sb_records(
    *,
    product_atoms: tuple[Any, ...],
    target_identities: tuple[tuple[str, str, int | None, str, str], ...],
    source_protein_pdb: Path,
    settings: InterchangeParameterizationSettings | None,
) -> tuple[AtomPartialChargeRecord, ...]:
    """Extract ff14SB-style charges from the prepared source protein."""
    source_charges = _source_protein_charges(source_protein_pdb, settings=settings)
    source_by_residue_atom = {
        (chain_id, residue_number, insertion_code, atom_name): charge
        for (
            chain_id,
            _residue_name,
            residue_number,
            insertion_code,
            atom_name,
        ), charge in source_charges.items()
    }
    product_atom_keys = {
        (
            atom.chain_id.strip(),
            atom.residue_number,
            atom.insertion_code.strip(),
            atom.atom_name.strip(),
        )
        for atom in product_atoms
        if atom.chain_id.strip() == "A"
    }
    records = []
    for identity in target_identities:
        chain_id, residue_name, residue_number, insertion_code, atom_name = identity
        if (
            chain_id != "A"
            or (chain_id, residue_number, insertion_code, atom_name) not in product_atom_keys
        ):
            continue
        charge = source_by_residue_atom.get((chain_id, residue_number, insertion_code, atom_name))
        if charge is None:
            continue
        records.append(
            AtomPartialChargeRecord(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                atom_name=atom_name,
                charge_e=charge,
                source="production:ff14SB-prepared-source-protein",
                source_role="protein_ff14sb",
            )
        )
    return tuple(records)


def _source_protein_charges(
    source_protein_pdb: Path,
    *,
    settings: InterchangeParameterizationSettings | None,
) -> dict[tuple[str, str, int | None, str, str], float]:
    """Parameterize the source protein and return charges by residue atom identity."""
    from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor

    ingested = PabloIngestor(policy=None).ingest_structure(source_protein_pdb)
    if not ingested.success or ingested.topology is None:
        raise ValueError(
            f"Could not ingest source protein for ff14SB charges: {source_protein_pdb}"
        )
    force_field = load_combined_smirnoff_force_field(settings)
    interchange = force_field.create_interchange(ingested.topology)
    charges = _interchange_charges_by_atom_index(interchange)
    atoms = tuple(_topology_atoms(ingested.topology))
    if len(charges) != len(atoms):
        raise ValueError(
            "Could not extract a complete ff14SB charge vector from source protein Interchange: "
            f"{len(charges)} charges for {len(atoms)} atom(s)"
        )
    return {_atom_identity(atom): charges[index] for index, atom in enumerate(atoms)}


def _polymer_template_records(specs: Sequence[Any]) -> tuple[AtomPartialChargeRecord, ...]:
    """Transfer attached-polymer charges from existing charged source templates."""
    records: list[AtomPartialChargeRecord] = []
    for spec in specs:
        sidecar = _source_sdf_path(spec)
        if sidecar is None:
            continue
        generated_fragment = getattr(spec, "generated_fragment", None)
        if generated_fragment is None:
            raise ValueError(
                "Attached polymer charge transfer requires generated_fragment metadata"
            )
        template_charges = _charged_sdf_atom_charges(sidecar)
        fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
        if len(template_charges) != len(fragment_atoms):
            raise ValueError(
                "Attached polymer charged SDF atom count does not match generated fragment: "
                f"{len(template_charges)} charges vs {len(fragment_atoms)} atom(s) for {sidecar}"
            )
        mappings = getattr(spec, "product_residue_mappings", {}) or {}
        leaving_names = {
            str(name).strip() for name in getattr(generated_fragment, "leaving_atom_names", ())
        }
        for atom, charge in zip(fragment_atoms, template_charges, strict=True):
            atom_name = str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")).strip()
            if atom_name in leaving_names:
                continue
            residue_number = int(getattr(atom, "residue_number", 0))
            mapped = _mapped_polymer_residue(atom, mappings)
            records.append(
                AtomPartialChargeRecord(
                    chain_id=str(mapped.get("chain_id", "C")),
                    residue_name=str(mapped.get("residue_name", getattr(atom, "residue_name", ""))),
                    residue_number=int(mapped.get("residue_number", residue_number)),
                    insertion_code=str(mapped.get("insertion_code", "")),
                    atom_name=atom_name,
                    charge_e=charge,
                    source=f"production:attached-polymer-template:{sidecar}",
                    source_role="polymer_template",
                )
            )
    return tuple(records)


def _local_nagl_patch_records(
    specs: Sequence[Any],
) -> tuple[tuple[AtomPartialChargeRecord, ...], str]:
    """Generate local NHS-Lys product-state NAGL patch charges."""
    charged_patch, model_name = _charge_patch_molecule()
    patch_charges = _mapped_patch_charges(charged_patch)
    records: list[AtomPartialChargeRecord] = []
    for spec in specs:
        plan = getattr(spec, "resolved_plan", None)
        if plan is None:
            raise ValueError("NHS-Lys local charge patch requires a resolved attachment plan")
        protein_atom = getattr(plan, "protein_link_atom", None)
        modifier_atom = getattr(plan, "modifier_link_atom", None)
        if protein_atom is None or modifier_atom is None:
            raise ValueError("NHS-Lys local charge patch requires protein and modifier link atoms")
        records.append(
            _record_from_pdb_atom(
                protein_atom,
                charge=patch_charges[_PATCH_NZ_MAP],
                residue_name=str(getattr(plan, "protein_product_residue_name", "LYX") or "LYX"),
                source=f"{_BRIDGE_SOURCE}:{model_name}:patch-map-{_PATCH_NZ_MAP}",
            )
        )
        records.append(
            _record_from_pdb_atom(
                modifier_atom,
                charge=patch_charges[_PATCH_CARBONYL_C_MAP],
                residue_name=str(getattr(plan, "modifier_product_residue_name", "") or ""),
                source=f"{_BRIDGE_SOURCE}:{model_name}:patch-map-{_PATCH_CARBONYL_C_MAP}",
            )
        )
        oxygen = _modifier_carbonyl_oxygen(spec)
        if oxygen is not None:
            records.append(
                _record_from_pdb_atom(
                    oxygen,
                    charge=patch_charges[_PATCH_CARBONYL_O_MAP],
                    residue_name=str(getattr(plan, "modifier_product_residue_name", "") or ""),
                    source=f"{_BRIDGE_SOURCE}:{model_name}:patch-map-{_PATCH_CARBONYL_O_MAP}",
                )
            )
    return tuple(records), model_name


def _charge_patch_molecule() -> tuple[Any, str]:
    """Charge the local product-state amide patch with NAGL."""
    from openff.toolkit import Molecule
    from rdkit import Chem

    from polyzymd.utils.charging import NAGLCharger

    rdmol = Chem.AddHs(Chem.MolFromSmiles(_PATCH_MAPPED_SMILES))
    molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for rd_atom, off_atom in zip(rdmol.GetAtoms(), molecule.atoms, strict=True):
        map_number = rd_atom.GetAtomMapNum()
        if map_number:
            off_atom.metadata["atom_map"] = map_number
    model_name = os.environ.get(PATCH_NAGL_MODEL_ENV, DEFAULT_PATCH_NAGL_MODEL)
    charger = NAGLCharger(model_name=model_name)
    charged = charger.charge_molecule(molecule)
    if getattr(charged, "partial_charges", None) is None:
        raise ValueError(f"NAGL model {model_name} did not return local patch partial charges")
    return charged, model_name


def _mapped_patch_charges(charged_patch: Any) -> dict[int, float]:
    """Return local patch charges by atom-map number."""
    charges = _charge_values(getattr(charged_patch, "partial_charges", ()))
    mapped: dict[int, float] = {}
    for atom, charge in zip(tuple(charged_patch.atoms), charges, strict=True):
        map_number = (getattr(atom, "metadata", None) or {}).get("atom_map")
        if map_number:
            mapped[int(map_number)] = charge
    required = {_PATCH_NZ_MAP, _PATCH_CARBONYL_C_MAP, _PATCH_CARBONYL_O_MAP}
    if not required.issubset(mapped):
        raise ValueError(f"NAGL patch molecule is missing mapped charges for {sorted(required)}")
    return mapped


def _modifier_carbonyl_oxygen(spec: Any) -> Any | None:
    """Return the first oxygen atom in the modifier link residue, if available."""
    plan = getattr(spec, "resolved_plan", None)
    modifier_atom = getattr(plan, "modifier_link_atom", None)
    fragment = getattr(spec, "generated_fragment", None)
    if modifier_atom is None or fragment is None:
        return None
    residue_number = getattr(modifier_atom, "residue_number", None)
    for atom in tuple(getattr(fragment, "atoms", ()) or ()):
        if getattr(atom, "residue_number", None) != residue_number:
            continue
        element = str(getattr(atom, "element", "") or "").strip().upper()
        if element == "O":
            return atom
    return None


def _record_from_pdb_atom(
    atom: Any,
    *,
    charge: float,
    residue_name: str,
    source: str,
) -> AtomPartialChargeRecord:
    """Build one charge record from a PDB-like atom."""
    return AtomPartialChargeRecord(
        chain_id=str(getattr(atom, "chain_id", "") or ""),
        residue_name=residue_name or str(getattr(atom, "residue_name", "") or ""),
        residue_number=getattr(atom, "residue_number", None),
        insertion_code=str(getattr(atom, "insertion_code", "") or ""),
        atom_name=str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")),
        charge_e=charge,
        source=source,
        source_role="local_nagl_patch",
    )


def _source_sdf_path(spec: Any) -> Path | None:
    """Return the production charged SDF path from an attachment spec."""
    sidecars = getattr(spec, "source_sidecars", {}) or {}
    if not isinstance(sidecars, Mapping):
        return None
    path = sidecars.get("charged_sdf")
    if path is not None:
        return Path(path)
    raw_path = sidecars.get("bond_sdf") or sidecars.get("sdf")
    if raw_path is not None:
        raise ValueError(
            "Attached polymer production charge transfer requires source_sidecars['charged_sdf']; "
            f"refusing raw bond/geometry SDF {raw_path} as a partial-charge source"
        )
    return None


def _charged_sdf_atom_charges(path: Path) -> tuple[float, ...]:
    """Read partial charges from a production charged SDF via PolyzyMD loaders."""
    from polyzymd.utils import get_largest_offmol, topology_from_sdf

    topology = topology_from_sdf(path)
    molecule = get_largest_offmol(topology)
    partial_charges = getattr(molecule, "partial_charges", None)
    if partial_charges is None:
        raise ValueError(
            "Attached polymer SDF does not contain partial charges; refusing to use AM1-BCC, "
            f"Gasteiger, or formal fallback for {path}"
        )
    return _charge_values(partial_charges)


def _mapped_polymer_residue(atom: Any, mappings: Mapping[str, Any]) -> dict[str, Any]:
    """Return mapped product residue identity for a generated fragment atom."""
    source_number = str(getattr(atom, "residue_number", "") or "")
    mapping = mappings.get(source_number) or mappings.get(
        f"{source_number}{getattr(atom, 'insertion_code', '')}"
    )
    if not mapping:
        return {
            "chain_id": getattr(atom, "chain_id", "C") or "C",
            "residue_name": getattr(atom, "residue_name", ""),
            "residue_number": getattr(atom, "residue_number", None),
            "insertion_code": getattr(atom, "insertion_code", "") or "",
        }
    return {
        "chain_id": mapping.get("target_chain", "C"),
        "residue_name": mapping.get("target_residue_name", getattr(atom, "residue_name", "")),
        "residue_number": mapping.get(
            "target_residue_number", getattr(atom, "residue_number", None)
        ),
        "insertion_code": mapping.get("target_insertion_code", ""),
    }


def _merge_records(
    records: dict[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    new_records: tuple[AtomPartialChargeRecord, ...],
    *,
    replace: bool = False,
) -> None:
    """Merge charge records with duplicate protection."""
    for record in new_records:
        key = record.identity_key
        if key in records and not replace:
            raise ValueError(f"Duplicate product-state charge source for {key}")
        records[key] = record


def _normalization_target(
    records: Mapping[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
) -> tuple[str, str, int | None, str, str] | None:
    """Return a conservative attached-polymer atom for total-charge closure."""
    for key, record in reversed(tuple(records.items())):
        if record.source_role == "polymer_template":
            return key
    return None


def _topology_atoms(topology: Any) -> tuple[Any, ...]:
    """Return atoms from an OpenFF topology in topology order."""
    atoms = getattr(topology, "atoms", None)
    if atoms is not None:
        return tuple(atoms)
    flattened = []
    for molecule in tuple(getattr(topology, "molecules", ()) or ()):
        flattened.extend(tuple(getattr(molecule, "atoms", ()) or ()))
    return tuple(flattened)


def _interchange_charges_by_atom_index(interchange: Any) -> dict[int, float]:
    """Extract electrostatic charges from an Interchange object."""
    collection = interchange["Electrostatics"]
    charges = getattr(collection, "charges", None)
    if charges is None:
        potentials = getattr(collection, "potentials", {})
        charges = {
            key: getattr(potential, "parameters", {}).get("charge")
            for key, potential in potentials.items()
        }
    values: dict[int, float] = {}
    for key, charge in dict(charges).items():
        atom_indices = tuple(getattr(key, "atom_indices", ()) or ())
        if len(atom_indices) != 1:
            continue
        values[int(atom_indices[0])] = _quantity_to_float(charge)
    return values


def _atom_identity(atom: Any) -> tuple[str, str, int | None, str, str]:
    """Return chain/residue/atom identity for an OpenFF-like atom."""
    metadata = getattr(atom, "metadata", None) or {}
    return (
        str(_metadata_value(metadata, "chain_id", "chain", default="") or "").strip(),
        str(_metadata_value(metadata, "residue_name", "residue", default="") or "").strip().upper(),
        _optional_int(_metadata_value(metadata, "residue_number", "residue_id", default=None)),
        str(_metadata_value(metadata, "insertion_code", default="") or "").strip(),
        str(
            _metadata_value(metadata, "atom_name", "name", default=getattr(atom, "name", "")) or ""
        ).strip(),
    )


def _metadata_value(metadata: Any, *names: str, default: Any = None) -> Any:
    """Return the first populated metadata value."""
    for name in names:
        value = metadata.get(name) if isinstance(metadata, dict) else getattr(metadata, name, None)
        if value not in (None, ""):
            return value
    return default


def _charge_values(quantity: Any) -> tuple[float, ...]:
    """Return elementary-charge floats from a charge quantity."""
    if hasattr(quantity, "m_as"):
        return tuple(float(value) for value in quantity.m_as("elementary_charge"))
    if hasattr(quantity, "value_in_unit"):
        from openff.units.openmm import to_openmm
        from openmm import unit

        try:
            values = quantity.value_in_unit(unit.elementary_charge)
        except Exception:  # noqa: BLE001 - support toolkit quantity variants
            values = to_openmm(quantity).value_in_unit(unit.elementary_charge)
        return tuple(float(value) for value in values)
    return tuple(float(value) for value in quantity)


def _quantity_to_float(quantity: Any) -> float:
    """Return one elementary-charge float from a quantity-like value."""
    if hasattr(quantity, "m_as"):
        return float(quantity.m_as("elementary_charge"))
    if hasattr(quantity, "value_in_unit"):
        from openff.units.openmm import to_openmm
        from openmm import unit

        try:
            return float(quantity.value_in_unit(unit.elementary_charge))
        except Exception:  # noqa: BLE001 - support toolkit quantity variants
            return float(to_openmm(quantity).value_in_unit(unit.elementary_charge))
    return float(quantity)


def _formal_charge_value(formal_charge: Any) -> float:
    """Return formal charge as an elementary-charge float."""
    return _quantity_to_float(formal_charge) if formal_charge is not None else 0.0


def _optional_int(value: Any) -> int | None:
    """Return an optional integer from metadata."""
    if value in (None, ""):
        return None
    return int(value)


def _format_identity(identity: tuple[str, str, int | None, str, str]) -> str:
    """Format an atom identity for diagnostics."""
    chain_id, residue_name, residue_number, insertion_code, atom_name = identity
    return (
        f"chain {chain_id or '?'} residue {residue_name or '?'} "
        f"{residue_number if residue_number is not None else '?'}{insertion_code} "
        f"atom {atom_name or '?'}"
    )
